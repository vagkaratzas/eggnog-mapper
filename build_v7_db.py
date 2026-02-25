#!/usr/bin/env python3
"""
build_v7_db.py — Convert eggNOG v7 flat files to v5-compatible eggnog.db

Usage:
    python3 build_v7_db.py \\
        --v7_dir  /path/to/v7/ \\
        --out_dir /path/to/v7_parsed/ \\
        --v5_dir  /path/to/v5/ \\
        --skip_events

The script produces three files in --out_dir:
    eggnog.db                    SQLite with version/og/prots/event tables
    eggnog.taxa.db               copied from v5 (pure NCBI taxonomy)
    eggnog.taxa.db.traverse.pkl  copied from v5
    eggnog_proteins.dmnd         built by diamond makedb (unless --skip_diamond)
"""

import argparse
import gzip
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ──────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert eggNOG v7 flat files to v5-compatible eggnog.db"
    )
    p.add_argument("--v7_dir", required=True,
                   help="Path to v7 flat files directory")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for parsed DB files")
    p.add_argument("--v5_dir",
                   help="Path to v5 dir (for taxa DB files to copy; "
                        "defaults to sibling 'v5' dir)")
    p.add_argument("--skip_diamond", action="store_true",
                   help="Skip diamond makedb step")
    p.add_argument("--skip_events", action="store_true",
                   help="Skip event table build (orthoindex stays NULL)")
    p.add_argument("--include_uncertain", action="store_true",
                   help="Include !-flagged OG assignments in prots.ogs")
    p.add_argument("--batch_size", type=int, default=50000,
                   help="SQLite executemany batch size (default: 50000)")
    p.add_argument("--threads", type=int, default=8,
                   help="diamond makedb thread count (default: 8)")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing eggnog.db")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# Field parsers
# ──────────────────────────────────────────────────────────────────────────────

def parse_kegg_ko(raw):
    """'K15082|7.14;K02315|8.33' → 'ko:K15082,ko:K02315'"""
    if not raw:
        return None
    parts = [e.split('|')[0] for e in raw.split(';') if e.strip()]
    return ','.join('ko:' + p for p in parts) if parts else None


def parse_pname(raw):
    """'RAD7|7.14;dnaI|16.67' → 'RAD7' (first/highest-score entry)"""
    if not raw:
        return None
    first = raw.split(';')[0].split('|')[0].strip()
    return first or None


def parse_gos_v7(raw):
    """'GO:0007010|30.77;GO:0000278|2.86' → 'u|GO:0007010|IEA,u|GO:0000278|IEA'

    Reformats v7 GO data (GO:id|score) into the three-part format expected by
    eggnogmapper/annotation/annota.py::parse_gos() (gocat|GO:id|evidence).
    """
    if not raw:
        return None
    parts = ['u|' + e.split('|')[0] + '|IEA' for e in raw.split(';') if e.strip()]
    return ','.join(parts) if parts else None


def split_og_assignment(og_str):
    """Parse 'ClusterName@taxid|ShortCode[!*]' into components.

    Returns (cluster_name, taxid, is_uncertain).
    Returns (None, None, False) on parse error.

    The '!' suffix marks uncertain OG assignments; '*' marks the seed/rep OG
    (not uncertain).  Both are stripped before further parsing.
    """
    is_uncertain = og_str.endswith('!')
    clean = og_str.rstrip('*!')
    at_idx = clean.find('@')
    if at_idx < 0:
        return None, None, False
    cluster_name = clean[:at_idx]
    rest = clean[at_idx + 1:]   # 'taxid|ShortCode'
    pipe_idx = rest.find('|')
    taxid = rest[:pipe_idx] if pipe_idx >= 0 else rest
    return cluster_name, taxid, is_uncertain


# ──────────────────────────────────────────────────────────────────────────────
# NHX tree parser (used only in Phase 3)
# ──────────────────────────────────────────────────────────────────────────────

class _NHXNode:
    """Lightweight Newick/NHX tree node."""
    __slots__ = ('children', 'name', 'nhx', '_leaves')

    def __init__(self):
        self.children = []
        self.name = ''
        self.nhx = {}
        self._leaves = None

    def leaves(self):
        """Return frozenset of leaf names (cached)."""
        if self._leaves is not None:
            return self._leaves
        if not self.children:
            self._leaves = frozenset((self.name,)) if self.name else frozenset()
        else:
            s = set()
            for c in self.children:
                s.update(c.leaves())
            self._leaves = frozenset(s)
        return self._leaves


def _parse_nhx_tree(s):
    """Parse a Newick/NHX string and return the root _NHXNode.

    Handles the NHX extension format used by eggNOG v7:
        node_label:branch_length[&&NHX:key=val:key=val]
    where internal nodes also carry a support value as label.
    """
    pos = [0]
    n = len(s)

    def skip_ws():
        while pos[0] < n and s[pos[0]] in ' \t\n\r':
            pos[0] += 1

    def read_nhx_block():
        """Consume '[&&NHX:...]', return dict of parsed key=value pairs."""
        pos[0] += 1  # skip '['
        depth, buf = 1, []
        while pos[0] < n:
            c = s[pos[0]]; pos[0] += 1
            if c == '[':
                depth += 1; buf.append(c)
            elif c == ']':
                depth -= 1
                if depth == 0:
                    break
                buf.append(c)
            else:
                buf.append(c)
        content = ''.join(buf)
        nhx = {}
        if content.startswith('&&NHX:'):
            for kv in content[6:].split(':'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    nhx[k] = v
        return nhx

    def read_label():
        """Read a node label, stopping at Newick structural characters or '['."""
        buf = []
        while pos[0] < n and s[pos[0]] not in ',:()[];[':
            buf.append(s[pos[0]]); pos[0] += 1
        return ''.join(buf)

    def read_branch():
        """Skip past a branch length (after ':'), stopping before '[' or Newick delimiters."""
        pos[0] += 1  # skip ':'
        while pos[0] < n and s[pos[0]] not in ',)[];[':
            pos[0] += 1

    def read_subtree():
        skip_ws()
        node = _NHXNode()
        if pos[0] < n and s[pos[0]] == '(':
            pos[0] += 1  # consume '('
            while True:
                skip_ws()
                node.children.append(read_subtree())
                skip_ws()
                if pos[0] < n and s[pos[0]] == ',':
                    pos[0] += 1
                elif pos[0] < n and s[pos[0]] == ')':
                    pos[0] += 1
                    break
                else:
                    break  # malformed tree; stop gracefully
            node.name = read_label()
            skip_ws()
            if pos[0] < n and s[pos[0]] == ':':
                read_branch()
            skip_ws()
            if pos[0] < n and s[pos[0]] == '[':
                node.nhx = read_nhx_block()
        else:
            node.name = read_label()
            skip_ws()
            if pos[0] < n and s[pos[0]] == ':':
                read_branch()
            skip_ws()
            if pos[0] < n and s[pos[0]] == '[':
                node.nhx = read_nhx_block()
        return node

    return read_subtree()


def collect_speciation_events(root, og_name):
    """Iteratively walk the tree and yield speciation events.

    Yields (lca_taxid, og_name, side1_csv, side2_csv) for every internal node
    annotated with e=S (speciation).  Uses an explicit stack to avoid Python
    recursion limits on deep trees.
    """
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            stack.append(child)
        if node.nhx.get('e') == 'S' and len(node.children) >= 2:
            lca = node.nhx.get('lca', '')
            side1 = ','.join(sorted(node.children[0].leaves()))
            side2 = ','.join(sorted(node.children[1].leaves()))
            yield lca, og_name, side1, side2


# ──────────────────────────────────────────────────────────────────────────────
# Phase 0 — Setup
# ──────────────────────────────────────────────────────────────────────────────

_SCHEMA = """
CREATE TABLE IF NOT EXISTS version (
    version VARCHAR(16) PRIMARY KEY
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS og (
    og               VARCHAR(64),
    level            VARCHAR(16),
    nm               INTEGER,
    description      TEXT,
    COG_categories   VARCHAR(8),
    PRIMARY KEY (og, level)
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS prots (
    name             VARCHAR(64) PRIMARY KEY,
    bigg_reaction    VARCHAR(32),
    gos              TEXT,
    pfam             TEXT,
    pname            VARCHAR(64),
    ogs              TEXT,
    orthoindex       TEXT,
    kegg_ko          VARCHAR(512),
    kegg_cog         VARCHAR(256),
    kegg_disease     VARCHAR(256),
    kegg_ec          VARCHAR(256),
    kegg_brite       VARCHAR(256),
    kegg_rclass      VARCHAR(256),
    kegg_tc          VARCHAR(256),
    kegg_cazy        VARCHAR(256),
    kegg_pathway     VARCHAR(256),
    kegg_module      VARCHAR(256),
    kegg_reaction    VARCHAR(256),
    kegg_go          VARCHAR(256),
    kegg_drug        VARCHAR(256),
    kegg_pubmed      TEXT,
    kegg_network     TEXT
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS event (
    i       INTEGER PRIMARY KEY,
    level   VARCHAR(16),
    og      VARCHAR(64),
    side1   TEXT,
    side2   TEXT
);
"""


def phase0_setup(args, out_dir, v5_dir):
    """Create output directory, copy taxa files, create DB schema."""
    out_dir.mkdir(parents=True, exist_ok=True)

    db_path = out_dir / 'eggnog.db'
    if db_path.exists():
        if args.force:
            db_path.unlink()
            print(f"Removed existing {db_path}")
        else:
            print(f"ERROR: {db_path} already exists. Use --force to overwrite.",
                  file=sys.stderr)
            sys.exit(1)

    # Copy taxonomy files from v5 (pure NCBI taxonomy, version-independent)
    if v5_dir:
        for fname in ('eggnog.taxa.db', 'eggnog.taxa.db.traverse.pkl'):
            src = v5_dir / fname
            dst = out_dir / fname
            if src.exists():
                if not dst.exists():
                    shutil.copy2(src, dst)
                    print(f"Copied {fname}")
                else:
                    print(f"Skipping {fname} (already present in out_dir)")
            else:
                print(f"WARNING: {src} not found, skipping", file=sys.stderr)
    else:
        print("WARNING: --v5_dir not provided; taxa DB files not copied.",
              file=sys.stderr)

    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=100000")
    conn.executescript(_SCHEMA)
    conn.execute("INSERT OR IGNORE INTO version (version) VALUES ('7.0')")
    conn.commit()
    print(f"Created schema in {db_path}")
    return conn


# ──────────────────────────────────────────────────────────────────────────────
# Phase 1 — og table + annotation cache
# ──────────────────────────────────────────────────────────────────────────────

def phase1_og_table(args, conn, v7_dir, level_depth):
    """Stream e7.og_info_kegg_go.tsv.gz, populate og table, build og_cache.

    og_cache maps (cluster_name, taxid) → (kegg_ko, pname, gos) for all
    entries regardless of whether taxid is in LEVEL_DEPTH.  Only entries whose
    taxid IS in LEVEL_DEPTH are inserted into the og table.

    Returns og_cache dict (~400 MB RAM for full v7 dataset).
    """
    og_info_file = v7_dir / 'e7.og_info_kegg_go.tsv.gz'
    print(f"\n[Phase 1] Building og table from {og_info_file.name}")

    og_cache = {}
    og_batch = []
    n_inserted = 0
    n_skipped_level = 0

    with gzip.open(og_info_file, 'rt', encoding='utf-8', errors='replace') as f:
        for lineno, line in enumerate(f, 1):
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue

            # cols[1] = cluster_name, cols[2] = taxid — parsed directly
            cluster_name = cols[1].strip()
            taxid = cols[2].strip()
            if not cluster_name or not taxid:
                continue

            nm_raw = cols[3].strip()
            nm = int(nm_raw) if nm_raw.isdigit() else 0

            kegg_ko = parse_kegg_ko(cols[6].strip()) if len(cols) > 6 else None
            pname   = parse_pname(cols[7].strip())   if len(cols) > 7 else None
            gos     = parse_gos_v7(cols[8].strip())  if len(cols) > 8 else None

            og_cache[(cluster_name, taxid)] = (kegg_ko, pname, gos)

            if taxid in level_depth:
                og_batch.append((cluster_name, taxid, nm, '', ''))
                n_inserted += 1
            else:
                n_skipped_level += 1

            if len(og_batch) >= args.batch_size:
                conn.executemany(
                    "INSERT OR IGNORE INTO og "
                    "(og, level, nm, description, COG_categories) "
                    "VALUES (?,?,?,?,?)",
                    og_batch
                )
                conn.commit()
                og_batch.clear()
                print(f"  {n_inserted:,} og rows inserted (line {lineno:,})...",
                      end='\r')

    if og_batch:
        conn.executemany(
            "INSERT OR IGNORE INTO og "
            "(og, level, nm, description, COG_categories) "
            "VALUES (?,?,?,?,?)",
            og_batch
        )
        conn.commit()

    print(f"\n  og table : {n_inserted:,} rows inserted, "
          f"{n_skipped_level:,} skipped (taxid not in LEVEL_DEPTH)")
    print(f"  og_cache : {len(og_cache):,} entries")
    return og_cache


# ──────────────────────────────────────────────────────────────────────────────
# Phase 2 — prots table
# ──────────────────────────────────────────────────────────────────────────────

_PROTS_INSERT = (
    "INSERT OR IGNORE INTO prots "
    "(name, bigg_reaction, gos, pfam, pname, ogs, orthoindex, kegg_ko, "
    "kegg_cog, kegg_disease, kegg_ec, kegg_brite, kegg_rclass, kegg_tc, "
    "kegg_cazy, kegg_pathway, kegg_module, kegg_reaction, kegg_go, "
    "kegg_drug, kegg_pubmed, kegg_network) "
    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
)


def phase2_prots_table(args, conn, v7_dir, og_cache, level_depth):
    """Stream e7.protein_families.tsv.gz and populate prots table.

    For each cluster row:
      col[4]  comma-separated protein IDs → one prots row each
      col[6]  comma-separated OG assignments (ClusterName@taxid|ShortCode[!*])

    OG assignments flagged with '!' are skipped unless --include_uncertain.
    Only OG assignments whose taxid is in LEVEL_DEPTH contribute to ogs.
    Annotations (kegg_ko, pname, gos) are merged across all valid OG assignments.
    """
    pf_file = v7_dir / 'e7.protein_families.tsv.gz'
    print(f"\n[Phase 2] Building prots table from {pf_file.name}")

    prots_batch = []
    n_total = 0
    n_no_ogs = 0

    with gzip.open(pf_file, 'rt', encoding='utf-8', errors='replace') as f:
        for lineno, line in enumerate(f, 1):
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                continue

            protein_ids_raw = cols[4].strip()
            if not protein_ids_raw:
                continue
            protein_ids = [p.strip() for p in protein_ids_raw.split(',') if p.strip()]

            og_assignments_raw = cols[6].strip() if len(cols) > 6 else ''

            valid_ogs = []
            all_ko = []
            all_gos = []
            pname = None

            if og_assignments_raw:
                for og_str in og_assignments_raw.split(','):
                    og_str = og_str.strip()
                    if not og_str:
                        continue
                    cn, taxid, is_uncertain = split_og_assignment(og_str)
                    if cn is None:
                        continue
                    if is_uncertain and not args.include_uncertain:
                        continue

                    ann = og_cache.get((cn, taxid), (None, None, None))
                    ko, pn, go = ann
                    if ko:
                        all_ko.extend(ko.split(','))
                    if go:
                        all_gos.extend(go.split(','))
                    if pn and pname is None:
                        pname = pn

                    if taxid in level_depth:
                        valid_ogs.append(f"{cn}@{taxid}")

            if not valid_ogs:
                # No valid LEVEL_DEPTH OGs → annotator would skip anyway; skip to save space
                n_no_ogs += len(protein_ids)
                continue

            ogs_str = ','.join(valid_ogs)
            ko_str  = ','.join(dict.fromkeys(all_ko)) or None
            gos_str = ','.join(dict.fromkeys(all_gos)) or None

            for prot_id in protein_ids:
                prots_batch.append((
                    prot_id,   # name
                    None,      # bigg_reaction
                    gos_str,   # gos
                    None,      # pfam
                    pname,     # pname
                    ogs_str,   # ogs
                    None,      # orthoindex
                    ko_str,    # kegg_ko
                    None, None, None, None, None, None,  # kegg_cog…kegg_tc
                    None, None, None, None, None, None,  # kegg_cazy…kegg_drug
                    None, None,                          # kegg_pubmed, kegg_network
                ))
                n_total += 1

            if len(prots_batch) >= args.batch_size:
                conn.executemany(_PROTS_INSERT, prots_batch)
                conn.commit()
                prots_batch.clear()
                print(f"  {n_total:,} proteins inserted (line {lineno:,})...",
                      end='\r')

    if prots_batch:
        conn.executemany(_PROTS_INSERT, prots_batch)
        conn.commit()

    print(f"\n  prots table  : {n_total:,} proteins inserted")
    print(f"  No valid OGs : {n_no_ogs:,} proteins")
    return n_total


# ──────────────────────────────────────────────────────────────────────────────
# Phase 3 — event table + orthoindex
# ──────────────────────────────────────────────────────────────────────────────

def phase3_events(args, conn, v7_dir):
    """Parse NHX trees from e7.trees.tsv.gz, populate event table, update orthoindex.

    e7.trees.tsv.gz columns (0-indexed):
        [0] cluster_id   [3] NHX Newick tree string

    For each e=S (speciation) internal node:
        i       = global event counter
        level   = lca taxid from [&&NHX:lca=TAXID]
        og      = cluster_id
        side1   = comma-sep sorted leaf names from child[0]
        side2   = comma-sep sorted leaf names from child[1]

    After all trees: orthoindex in prots is set to a comma-separated list of
    event IDs for all speciation events in which the protein appears as a leaf.
    """
    trees_file = v7_dir / 'e7.trees.tsv.gz'
    print(f"\n[Phase 3] Building event table from {trees_file.name}")

    event_batch = []
    event_counter = 0
    protein_to_events = {}   # protein_id → [i, ...]
    n_events = 0
    n_parse_errors = 0

    with gzip.open(trees_file, 'rt', encoding='utf-8', errors='replace') as f:
        for lineno, line in enumerate(f, 1):
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue
            og_name = cols[0].strip()
            tree_str = cols[3].strip()
            if not tree_str:
                continue

            try:
                root = _parse_nhx_tree(tree_str)
            except Exception as exc:
                n_parse_errors += 1
                if n_parse_errors <= 10:
                    print(f"\n  WARNING: parse error line {lineno} "
                          f"(og={og_name}): {exc}", file=sys.stderr)
                continue

            for lca, og, side1, side2 in collect_speciation_events(root, og_name):
                i = event_counter
                event_counter += 1
                event_batch.append((i, lca, og, side1, side2))
                n_events += 1

                for prot in side1.split(','):
                    if prot:
                        protein_to_events.setdefault(prot, []).append(i)
                for prot in side2.split(','):
                    if prot:
                        protein_to_events.setdefault(prot, []).append(i)

            if len(event_batch) >= args.batch_size:
                conn.executemany(
                    "INSERT OR IGNORE INTO event (i, level, og, side1, side2) "
                    "VALUES (?,?,?,?,?)",
                    event_batch
                )
                conn.commit()
                event_batch.clear()
                print(f"  {n_events:,} events (line {lineno:,})...", end='\r')

    if event_batch:
        conn.executemany(
            "INSERT OR IGNORE INTO event (i, level, og, side1, side2) "
            "VALUES (?,?,?,?,?)",
            event_batch
        )
        conn.commit()

    if n_parse_errors:
        print(f"\n  WARNING: {n_parse_errors:,} trees failed to parse",
              file=sys.stderr)
    print(f"\n  event table : {n_events:,} events")
    print(f"  Updating orthoindex for {len(protein_to_events):,} proteins...")

    update_batch = []
    for prot_id, event_ids in protein_to_events.items():
        update_batch.append((','.join(str(x) for x in event_ids), prot_id))
        if len(update_batch) >= args.batch_size:
            conn.executemany(
                "UPDATE prots SET orthoindex=? WHERE name=?",
                update_batch
            )
            conn.commit()
            update_batch.clear()

    if update_batch:
        conn.executemany(
            "UPDATE prots SET orthoindex=? WHERE name=?",
            update_batch
        )
        conn.commit()

    print("  orthoindex updated.")
    return n_events


# ──────────────────────────────────────────────────────────────────────────────
# Phase 4 — diamond makedb
# ──────────────────────────────────────────────────────────────────────────────

def phase4_diamond(args, v7_dir, out_dir):
    """Build DIAMOND database from e7.proteins.fa.gz."""
    fasta = v7_dir / 'e7.proteins.fa.gz'
    dmnd  = out_dir / 'eggnog_proteins.dmnd'
    print(f"\n[Phase 4] Running diamond makedb → {dmnd.name}")
    result = subprocess.run(
        ['diamond', 'makedb',
         '--in', str(fasta),
         '--db', str(dmnd),
         '--threads', str(args.threads)],
        check=False
    )
    if result.returncode != 0:
        print("WARNING: diamond makedb exited with non-zero status.",
              file=sys.stderr)
    else:
        print("  diamond makedb complete.")


# ──────────────────────────────────────────────────────────────────────────────
# Phase 5 — finalize
# ──────────────────────────────────────────────────────────────────────────────

def phase5_finalize(conn):
    """Switch WAL → standard journal mode, print row counts, close."""
    print("\n[Phase 5] Finalizing SQLite database...")
    conn.execute("PRAGMA journal_mode=DELETE")
    for table in ('version', 'og', 'prots', 'event'):
        count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table:10s}: {count:,} rows")
    conn.close()
    print("  Done.")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # Import LEVEL_DEPTH from eggnog-mapper package in same directory
    repo_root = Path(__file__).parent
    sys.path.insert(0, str(repo_root))
    try:
        from eggnogmapper.annotation.tax_scopes.vars import LEVEL_DEPTH
    except ImportError as e:
        print(f"ERROR: cannot import LEVEL_DEPTH from eggnogmapper: {e}",
              file=sys.stderr)
        print("Run this script from the eggnog-mapper repository root.",
              file=sys.stderr)
        sys.exit(1)

    v7_dir  = Path(args.v7_dir)
    out_dir = Path(args.out_dir)
    v5_dir  = Path(args.v5_dir) if args.v5_dir else v7_dir.parent / 'v5'

    for d, label in ((v7_dir, '--v7_dir'),):
        if not d.is_dir():
            print(f"ERROR: {label} directory not found: {d}", file=sys.stderr)
            sys.exit(1)

    print(f"v7_dir  : {v7_dir}")
    print(f"out_dir : {out_dir}")
    print(f"v5_dir  : {v5_dir}")
    print(f"LEVEL_DEPTH contains {len(LEVEL_DEPTH):,} recognized taxids")

    # Phase 0 — setup
    conn = phase0_setup(args, out_dir, v5_dir)

    # Phase 1 — og table + annotation cache
    og_cache = phase1_og_table(args, conn, v7_dir, LEVEL_DEPTH)

    # Phase 2 — prots table
    phase2_prots_table(args, conn, v7_dir, og_cache, LEVEL_DEPTH)
    del og_cache  # free ~400 MB

    # Phase 3 — event table (optional)
    if not args.skip_events:
        phase3_events(args, conn, v7_dir)
    else:
        print("\n[Phase 3] Skipped (--skip_events).")

    # Phase 4 — diamond makedb (optional)
    if not args.skip_diamond:
        phase4_diamond(args, v7_dir, out_dir)
    else:
        print("\n[Phase 4] Skipped (--skip_diamond).")

    # Phase 5 — finalize
    phase5_finalize(conn)


if __name__ == '__main__':
    main()
