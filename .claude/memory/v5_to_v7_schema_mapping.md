# eggNOG v5 → v7 Database Schema Mapping

## Overview

eggNOG v5 ships pre-built SQLite (`eggnog.db`) and DIAMOND (`eggnog_proteins.dmnd`)
databases. eggNOG v7 ships raw flat files (TSV/FASTA) that must be ingested by
`build_v7_db.py` to produce the equivalent layout.

---

## `eggnog.db` tables

### `version`

| Column | Type | v5 value | v7 value |
|---|---|---|---|
| `version` | VARCHAR(16) PK | `'5.0'` | `'7.0'` |

No schema change.

---

### `og`

**v5 DDL:**
```sql
CREATE TABLE "og" (
    og             VARCHAR(16),
    level          VARCHAR(16),
    nm             INTEGER,
    description    TEXT,
    COG_categories VARCHAR(8),
    PRIMARY KEY (og, level)
) WITHOUT ROWID
```

**v7 DDL (build_v7_db.py):**
```sql
CREATE TABLE IF NOT EXISTS og (
    og             VARCHAR(64),   -- widened: v7 cluster names can exceed 16 chars
    level          VARCHAR(16),
    nm             INTEGER,
    description    TEXT,
    COG_categories VARCHAR(8),
    PRIMARY KEY (og, level)
) WITHOUT ROWID
```

**Note:** SQLite does not enforce VARCHAR lengths; this is cosmetic only.

**v7 source:** `e7.og_info_kegg_go.tsv.gz`

| DB column | v7 source column | Notes |
|---|---|---|
| `og` | col[1] — cluster_name | e.g. `YchF-GTPase_C` |
| `level` | col[2] — taxid string | e.g. `'131567'` |
| `nm` | col[3] — integer | cluster size |
| `description` | — | always `''`; v7 has no equivalent |
| `COG_categories` | — | always `''`; v7 has no equivalent |

**Filtering:** only rows where `taxid in LEVEL_DEPTH` (379-entry whitelist in
`eggnogmapper/annotation/tax_scopes/vars.py`) are inserted. ~14% of v7 taxids
qualify. Rows outside LEVEL_DEPTH are skipped to prevent `KeyError` crashes in
`parse_nogs()` (`tax_scopes.py:92`), which does a bare `LEVEL_DEPTH[nog_tax_id]`
lookup with no `.get()` fallback.

---

### `prots`

**v5 DDL:**
```sql
CREATE TABLE "prots" (
    name           VARCHAR(32) PRIMARY KEY,
    bigg_reaction  VARCHAR(32),
    gos            TEXT,
    pfam           TEXT,
    pname          VARCHAR(32),
    ogs            VARCHAR(256),
    orthoindex     VARCHAR(256),
    kegg_ko        VARCHAR(256),
    kegg_cog       VARCHAR(256),
    kegg_disease   VARCHAR(256),
    kegg_ec        VARCHAR(256),
    kegg_brite     VARCHAR(256),
    kegg_rclass    VARCHAR(256),
    kegg_tc        VARCHAR(256),
    kegg_cazy      VARCHAR(256),
    kegg_pathway   VARCHAR(256),
    kegg_module    VARCHAR(256),
    kegg_reaction  VARCHAR(256),
    kegg_go        VARCHAR(256),
    kegg_drug      VARCHAR(256),
    kegg_pubmed    TEXT,
    kegg_network   TEXT
) WITHOUT ROWID
```

**v7 DDL (build_v7_db.py) — changed columns highlighted:**
```sql
CREATE TABLE IF NOT EXISTS prots (
    name           VARCHAR(64) PRIMARY KEY,   -- widened (32→64)
    bigg_reaction  VARCHAR(32),
    gos            TEXT,
    pfam           TEXT,
    pname          VARCHAR(64),               -- widened (32→64)
    ogs            TEXT,                      -- widened (VARCHAR(256)→TEXT)
    orthoindex     TEXT,                      -- widened (VARCHAR(256)→TEXT)
    kegg_ko        VARCHAR(512),              -- widened (256→512)
    kegg_cog       VARCHAR(256),
    kegg_disease   VARCHAR(256),
    kegg_ec        VARCHAR(256),
    kegg_brite     VARCHAR(256),
    kegg_rclass    VARCHAR(256),
    kegg_tc        VARCHAR(256),
    kegg_cazy      VARCHAR(256),
    kegg_pathway   VARCHAR(256),
    kegg_module    VARCHAR(256),
    kegg_reaction  VARCHAR(256),
    kegg_go        VARCHAR(256),
    kegg_drug      VARCHAR(256),
    kegg_pubmed    TEXT,
    kegg_network   TEXT
) WITHOUT ROWID
```

**All width changes are cosmetic** — SQLite stores TEXT and VARCHAR(n) identically
(TEXT affinity). Any query or tool targeting the v5 schema works unchanged against
the v7-built DB.

**v7 source:** `e7.protein_families.tsv.gz` (protein→cluster mapping) joined
with annotation cache built from `e7.og_info_kegg_go.tsv.gz`

| DB column | v7 derivation | Notes |
|---|---|---|
| `name` | col[4] of protein_families — each protein ID | e.g. `199310.c0071` |
| `ogs` | col[6] of protein_families — OG assignments, filtered to LEVEL_DEPTH taxids, deduplicated | format: `ClusterName@taxid,...` |
| `pname` | col[7] of og_info (`RAD7\|7.14;...`) — first entry's name | requires freq≥2 across co-orthologs to appear in emapper output |
| `gos` | col[8] of og_info (`GO:0005975\|66.67;...`) — reformatted | stored as `u\|GO:id\|IEA`; filtered by emapper's `--go_evidence` (default excludes IEA; use `--go_evidence all`) |
| `kegg_ko` | col[6] of og_info (`K15082\|7.14;...`) — score stripped, `ko:` prefix added | stored as `ko:K15082,ko:K02315,...` |
| `kegg_ec` | — | NULL; v7 does not provide EC numbers per-protein |
| `kegg_pathway` | — | NULL |
| `kegg_module` | — | NULL |
| all other `kegg_*` | — | NULL |
| `pfam` | — | NULL |
| `bigg_reaction` | — | NULL |
| `orthoindex` | — | NULL unless `--skip_events` is omitted |

**Storage optimisation (hardcoded, no flag):** proteins whose every OG assignment
falls outside LEVEL_DEPTH are not inserted at all. Storing them with `ogs=NULL`
would be correct (emapper skips them) but wasteful — the full v7 has ~76M proteins
vs ~12M with the filter, saving ~60 GB.

---

### `event`

**v5 DDL:**
```sql
CREATE TABLE event (
    i       INTEGER PRIMARY KEY,
    level   VARCHAR(16),
    og      VARCHAR(16),
    side1   TEXT,
    side2   TEXT
)
```

**v7 DDL (build_v7_db.py):**
```sql
CREATE TABLE IF NOT EXISTS event (
    i       INTEGER PRIMARY KEY,
    level   VARCHAR(16),
    og      VARCHAR(64),   -- widened: v7 cluster names
    side1   TEXT,
    side2   TEXT
)
```

**v7 source:** `e7.trees.tsv.gz` — NHX Newick trees; `e=S` internal nodes are
speciation events. Built by `--skip_events` being absent. Skipped by default in
initial builds due to high time and disk cost.

When skipped, `prots.orthoindex` stays NULL and emapper falls back to seed-ortholog-
only for ortholog inference, which means `Preferred_name` is always empty in output.

---

## `eggnog.taxa.db` and `eggnog.taxa.db.traverse.pkl`

These are **not rebuilt from v7 flat files**. They are copied from the v5 directory.

- `eggnog.taxa.db` — full NCBI taxonomy (1,229,741 species rows). The v7 file
  `e7.taxid_info.tsv.gz` only covers the 22,893 taxa present in v7 eggNOG and lacks
  the `synonym` table and full tree needed for `get_descendant_taxa()`.
- `eggnog.taxa.db.traverse.pkl` — pre-computed pre-post-order traversal of the
  taxonomy tree; derived from `eggnog.taxa.db`, not from v7 files.

The NCBI taxonomy is version-independent between v5 and v7 for emapper's purposes.
The canonical way to build these from scratch is from NCBI's `taxdump.tar.gz`
(names.dmp, nodes.dmp, merged.dmp), not from eggNOG flat files.

---

## `eggnog_proteins.dmnd`

Built by `diamond makedb` from `e7.proteins.fa.gz`. No schema equivalent in v5;
this is a binary DIAMOND index. Contains all ~59.3M v7 proteins (full, unfiltered).

---

## Key format conversions (v7 → v5-compatible)

| Field | v7 raw format | v5-compatible format | Parser |
|---|---|---|---|
| GO terms | `GO:0005975\|66.67;GO:0000278\|2.86` | `u\|GO:0005975\|IEA,u\|GO:0000278\|IEA` | `parse_gos_v7()` |
| KEGG KO | `K15082\|7.14;K02315\|8.33` | `ko:K15082,ko:K02315` | `parse_kegg_ko()` |
| pname | `RAD7\|7.14;dnaI\|16.67` | `RAD7` (first entry) | `parse_pname()` |
| OG assignment | `ClusterName@taxid\|ShortCode[!*]` | `ClusterName@taxid` | `split_og_assignment()` |
