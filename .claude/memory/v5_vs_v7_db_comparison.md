# eggNOG v5 vs v7_parsed Database Comparison

Computed by querying both SQLite databases directly.
- v5 source: `/home/vangelis/Desktop/Databases/eggnogdb/v5/eggnog.db`
- v7 source: `/home/vangelis/Desktop/Databases/eggnogdb/v7_parsed/eggnog.db`
  (partial build — disk space exhausted mid-Phase 2; ~16% of full v7 proteins loaded)

---

## Row counts

| Table   | v5         | v7_parsed  |
|---------|------------|------------|
| `og`    | 4,415,548  | 352,360    |
| `prots` | 23,714,351 | 12,155,651 |
| `event` | 75,371,789 | 0          |

`event` is 0 in v7_parsed because the build was run with `--skip_events`.

---

## `og` table

### OG names
Completely non-overlapping naming systems — zero shared identifiers:
- v5 uses short alphanumeric codes: `1EQ0F`, `COG0001`, `KOG0001`, …
- v7 uses descriptive cluster names: `YchF-GTPase_C`, `DNA_pol_B_exo1`, …

### OG levels (taxids)

| | Count |
|---|---|
| v5 distinct levels | 379 |
| v7_parsed distinct levels | 282 |
| Shared (v7 ⊆ v5) | 282 |
| v5-only levels | 97 |
| v7-only levels | 0 |

All 282 v7_parsed levels are a strict subset of v5's 379. The 97 v5-only levels
fall into three groups:

**Virus clades** (dropped or reorganised in v7):
- root of viruses: Viruses (10239), Caudovirales (28883, 6,116 OGs),
  dsDNA viruses (35237, 6,987 OGs), Myoviridae (10662, 3,197 OGs),
  Siphoviridae (10699), Podoviridae (10744), Herpesvirales (548681),
  and ~15 smaller viral families

**Fine-grained eukaryote/prokaryote lineages** (not covered by v7_parsed LEVEL_DEPTH
filter — these 97 taxids are in v5's LEVEL_DEPTH but produced no OG rows in v7):
- Notable examples with high v5 OG counts:
  - root (taxid=1): 295,028 OGs — root is not in LEVEL_DEPTH whitelist
  - Insecta (50557): 22,816 OGs
  - Cercopithecoidea (314294): 20,011 OGs
  - Pseudonocardineae (85010): 21,205 OGs
  - Brassicales (3699): 17,310 OGs
  - Cytophagia (768503): 15,657 OGs
  - Testudines (8459): 16,870 OGs
  - Peronosporales (4776): 14,276 OGs
  - Drosophilidae (7214): 13,686 OGs
  - Nematocera (7148): 12,188 OGs

---

## `prots` table

| | Count | % of v5 | % of v7 |
|---|---|---|---|
| v5 total | 23,714,351 | — | — |
| v7_parsed total | 12,155,651 | — | — |
| **Shared** | **1,524,825** | **6.4%** | **12.5%** |
| v5-only | 22,189,526 | 93.6% | — |
| v7-only | 10,630,826 | — | 87.5% |

The low overlap (~1.5M proteins) is explained by three factors:

1. **v7_parsed is partial** — only ~16% of the full ~76M v7 proteins were inserted
   before disk space ran out (build stopped around line 1,174,313 of
   `e7.protein_families.tsv.gz`). The remaining ~64M v7 proteins were never loaded.

2. **v7 is a larger, newer database** — many v7 proteins represent organisms and
   assemblies added to NCBI since v5 was built, and therefore don't exist in v5.

3. **v5 covers different reference proteomes** — a large fraction of v5 proteins
   have identifiers that do not appear in v7, either because the sequences were
   removed, superseded, or renumbered between database versions.

A full v7 build (all ~76M proteins) would likely show a higher absolute overlap
with v5 but the percentage overlap would remain low given the scale difference.

---

## `event` table

| | v5 | v7_parsed |
|---|---|---|
| Rows | 75,371,789 | 0 |

v7_parsed was built with `--skip_events`. To populate this table, re-run
`build_v7_db.py` without `--skip_events` (requires substantial additional time
and ~50+ GB of extra disk space). Populating events enables:
- `prots.orthoindex` values (currently NULL for all v7 proteins)
- `Preferred_name` in emapper output (requires ≥2 co-orthologs)
- Proper one2one/one2many/many2many ortholog classification
