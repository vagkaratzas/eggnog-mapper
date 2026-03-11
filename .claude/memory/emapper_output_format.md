# eggNOG-mapper Output File Formats

emapper.py produces three output files per run, all tab-separated:

---

## 1. `.emapper.hits`

Raw DIAMOND search results. One row per hit (up to `--top N` hits per query).
Lines starting with `##` are comments; the header line starts with `#`.

**Header:** `#qseqid sseqid evalue bitscore qstart qend sstart send pident qcov scov`
(short format omits columns after bitscore; controlled by `--outfmt_short`)

| Column | Description |
|---|---|
| `qseqid` | Query sequence ID |
| `sseqid` | Subject (database hit) sequence ID — format `taxid.proteinid` |
| `evalue` | Alignment E-value |
| `bitscore` | Alignment bitscore |
| `qstart` | Query alignment start (1-based) |
| `qend` | Query alignment end |
| `sstart` | Subject alignment start |
| `send` | Subject alignment end |
| `pident` | Percent identity |
| `qcov` | Query coverage % |
| `scov` | Subject coverage % |

---

## 2. `.emapper.seed_orthologs`

One row per query, showing only the single best hit selected as the seed ortholog
after applying evalue/score/coverage filters.

**Header:** `#qseqid sseqid evalue bitscore qstart qend sstart send pident qcov scov`

Same columns as `.hits` but reduced to one row per query. This file is used as
input to `--annotate_hits_table` to re-run annotation without repeating the search.

---

## 3. `.emapper.annotations`

One row per query with functional annotations transferred from co-orthologs.
Lines starting with `##` are comments; the header line starts with `#`.

**Column order (21 columns):**

| # | Column | Source | Description |
|---|---|---|---|
| 1 | `query` | input | Query sequence ID |
| 2 | `seed_ortholog` | hits | Best database hit selected as seed (`taxid.proteinid`) |
| 3 | `evalue` | hits | E-value of the seed ortholog alignment |
| 4 | `score` | hits | Bitscore of the seed ortholog alignment |
| 5 | `eggNOG_OGs` | `prots.ogs` | All orthologous groups the seed belongs to, comma-separated, sorted narrowest→broadest. Format: `ClusterName@taxid\|TaxonName` |
| 6 | `max_annot_lvl` | best OG | Taxonomic level of the OG used as annotation basis (chosen by `--tax_scope` / `--tax_scope_mode`). Under v7 with no event table this falls back to `seed_ortholog@proteinid\|-` |
| 7 | `COG_category` | `og.COG_categories` | One-letter COG functional category codes (e.g. `J`, `KT`). Empty in v7 |
| 8 | `Description` | `og.description` | Free-text OG description. Empty in v7 |
| 9 | `Preferred_name` | `prots.pname` | Most common gene symbol among co-orthologs (requires ≥2 co-orthologs with same name; empty without event table) |
| 10 | `GOs` | `prots.gos` | GO term IDs, comma-separated (e.g. `GO:0006260,GO:0051604`). Filtered by `--go_evidence`; v7 stores only IEA so use `--go_evidence all` |
| 11 | `EC` | `prots.kegg_ec` | Enzyme Commission numbers (e.g. `2.7.7.7`). Empty in v7 |
| 12 | `KEGG_ko` | `prots.kegg_ko` | KEGG Orthology IDs, comma-separated (e.g. `ko:K02315,ko:K02324`) |
| 13 | `KEGG_Pathway` | `prots.kegg_pathway` | KEGG pathway map IDs (e.g. `map00010`). Empty in v7 |
| 14 | `KEGG_Module` | `prots.kegg_module` | KEGG module IDs (e.g. `M00001`). Empty in v7 |
| 15 | `KEGG_Reaction` | `prots.kegg_reaction` | KEGG reaction IDs. Empty in v7 |
| 16 | `KEGG_rclass` | `prots.kegg_rclass` | KEGG reaction class. Empty in v7 |
| 17 | `BRITE` | `prots.kegg_brite` | KEGG BRITE hierarchy entries. Empty in v7 |
| 18 | `KEGG_TC` | `prots.kegg_tc` | Transporter Classification IDs. Empty in v7 |
| 19 | `CAZy` | `prots.kegg_cazy` | Carbohydrate-Active enZyme family (e.g. `GH18`). Empty in v7 |
| 20 | `BiGG_Reaction` | `prots.bigg_reaction` | BiGG metabolic reaction IDs. Empty in v7 |
| 21 | `PFAMs` | `prots.pfam` | Pfam domain families (e.g. `PF00085`). Empty in v7 |

Missing values in any column are written as `-`.

An optional 22nd column `md5` is appended when `--md5` is passed.

---

## Annotation transfer logic

Annotations are **not** taken directly from the seed ortholog's own row in `prots`.
They are aggregated across all **co-orthologs** inferred via the event table:

1. `get_member_ogs()` → fetches `prots.ogs` for the seed ortholog
2. `parse_nogs()` → selects the best OG given `--tax_scope` / `--tax_scope_mode`
3. `get_member_orthologs()` → looks up speciation events (`prots.orthoindex` →
   `event` table) to find all co-orthologs in the best OG
4. If no events (empty event table or NULL orthoindex): falls back to `{seed}` only
5. `summarize_annotations()` → queries `prots` for all co-orthologs and aggregates:
   - `GOs`: union of GO IDs passing `--go_evidence` filter
   - `Preferred_name`: most common `pname` value, only if freq ≥ 2
   - All other fields: union of comma-separated values

---

## `--go_evidence` filter

Controls which GO evidence codes are included in the `GOs` column:

| Option | `target_go_ev` | `excluded_go_ev` | Effect |
|---|---|---|---|
| `non-electronic` (default) | None | `{ND, IEA}` | Excludes electronically inferred terms |
| `experimental` | `{EXP,IDA,IPI,IMP,IGI,IEP}` | `{ND, IEA}` | Only wet-lab evidence |
| `all` | None | None | All terms kept |

**v7 note:** all GO terms in the v7-built database are stored with evidence code
`IEA` (Inferred from Electronic Annotation). The default `non-electronic` mode
therefore produces no GO output. Use `--go_evidence all` with v7 data.

---

## Key `--tax_scope_mode` values

Controls which OG level is chosen as the annotation basis (`max_annot_lvl`):

| Mode | Behaviour |
|---|---|
| `broadest` | Broadest (shallowest) OG across all tax scope levels |
| `inner_broadest` | Broadest OG within the matching tax scope intersection |
| `inner_narrowest` (default) | Narrowest OG within the matching tax scope intersection |
| `narrowest` | Narrowest (deepest) OG regardless of tax scope |
