[![Build Status](https://travis-ci.com/eggnogdb/eggnog-mapper.svg?branch=master)](https://travis-ci.com/eggnogdb/eggnog-mapper)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=eggnog_mapper)

# Building eggNOG v7 databases (`build_v7_db.py`)

eggNOG v7 ships raw flat files (TSV/FASTA) rather than the pre-built SQLite and DIAMOND
databases that emapper expects. `build_v7_db.py` (repo root) ingests the v7 flat files and
produces a v5-compatible database layout that `emapper.py` can use directly via `--data_dir`.

## What it produces

| File | Description |
|---|---|
| `eggnog.db` | SQLite with `version`, `og`, `prots`, and `event` tables |
| `eggnog.taxa.db` | NCBI taxonomy DB (copied from an existing v5 directory) |
| `eggnog.taxa.db.traverse.pkl` | Taxonomy traversal cache (copied from v5) |
| `eggnog_proteins.dmnd` | DIAMOND database built from `e7.proteins.fa.gz` |

Required v7 source files (must be present in `--v7_dir`):

| File | Used for |
|---|---|
| `e7.og_info_kegg_go.tsv.gz` | `og` table + per-OG KEGG/GO annotation cache |
| `e7.protein_families.tsv.gz` | `prots` table (protein → OG mapping) |
| `e7.trees.tsv.gz` | `event` table (NHX speciation events; only with `--skip_events` omitted) |
| `e7.proteins.fa.gz` | DIAMOND database (skipped with `--skip_diamond`) |

## Usage

```bash
python3 build_v7_db.py \
    --v7_dir  /path/to/eggnog_v7/ \
    --out_dir /path/to/output/ \
    --v5_dir  /path/to/eggnog_v5/ \
    --skip_events          # recommended for initial builds; see caveats below
```

Then annotate with:

```bash
python3 emapper.py \
    -m diamond \
    -i proteins.fa \
    --data_dir /path/to/output/ \
    --go_evidence all \    # v7 GO terms are IEA-only; the default 'non-electronic' mode excludes them
    -o my_run
```

## Arguments

| Argument | Default | Description |
|---|---|---|
| `--v7_dir` | *(required)* | Directory containing eggNOG v7 flat files |
| `--out_dir` | *(required)* | Output directory for the parsed database files |
| `--v5_dir` | sibling `v5/` dir | Source of `eggnog.taxa.db` and `.traverse.pkl` (pure NCBI taxonomy, identical between versions) |
| `--skip_events` | off | Skip the `event` table build. Saves substantial time and disk I/O; see caveats |
| `--skip_diamond` | off | Skip `diamond makedb`; useful if the `.dmnd` file already exists |
| `--include_uncertain` | off | Include OG assignments flagged with `!` (uncertain) in `prots.ogs` |
| `--batch_size` | 50000 | SQLite `executemany` batch size; reduce if memory is tight |
| `--threads` | 8 | Thread count for `diamond makedb` |
| `--force` | off | Overwrite an existing `eggnog.db` in `--out_dir` |

## Caveats and known limitations

**Only OGs at LEVEL_DEPTH-recognized taxa are stored.** eggNOG-mapper's internal
`LEVEL_DEPTH` dictionary (379 entries) is used as a whitelist. The v7 database defines
OGs at ~13,000 unique taxids; only the ~14% that overlap LEVEL_DEPTH are written to
`prots.ogs`. Proteins whose every OG assignment falls outside LEVEL_DEPTH are silently
omitted from the `prots` table entirely.

**GO terms are IEA-only.** v7 reports GO term association percentages rather than
experimental evidence codes. The build script stores all GO terms with evidence code `IEA`.
emapper's default `--go_evidence non-electronic` filter excludes IEA terms; run with
`--go_evidence all` to retain them.

**`Preferred_name` requires the event table.** emapper populates the `Preferred_name`
output column only when the same name appears in ≥2 co-orthologs. Without the event table
(`--skip_events`), ortholog lookup falls back to the seed ortholog alone (1 protein), so
`Preferred_name` is always empty. Build without `--skip_events` to enable this field.

**COG categories and descriptions are empty.** v7 does not provide COG functional
category letters or free-text OG descriptions equivalent to v5. These columns are stored
as empty strings and will appear as `-` in emapper output.

**Disk space.** A full build (all proteins) requires ~120 GB: ~70 GB for `eggnog.db`
(without the event table) plus ~25 GB for `eggnog_proteins.dmnd`. The `event` table adds
further significant storage. Ensure sufficient free space before starting.

---

# Overview
**EggNOG-mapper** is a tool for fast functional annotation of novel sequences. It uses precomputed orthologous groups and phylogenies from the eggNOG database (http://eggnog5.embl.de) to transfer functional information from fine-grained orthologs only.

Common uses of eggNOG-mapper include the annotation of novel genomes, transcriptomes or even metagenomic gene catalogs.

The use of orthology predictions for functional annotation permits a higher precision than traditional homology searches (i.e. BLAST searches), as it avoids transferring annotations from close paralogs (duplicate genes with a higher chance of being involved in functional divergence).

Benchmarks comparing different eggNOG-mapper options against BLAST and InterProScan [can be found here](https://github.com/jhcepas/emapper-benchmark/blob/master/benchmark_analysis.ipynb).

EggNOG-mapper is also available as a public online resource: http://eggnog-mapper.embl.de

# Documentation
https://github.com/jhcepas/eggnog-mapper/wiki

# Citation

If you use this software, please cite:
```
[1] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
    prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
    Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
    Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293

[2] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
    orthology resource based on 5090 organisms and 2502 viruses. Jaime
    Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernández-Plaza, Sofia
    K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas Rattei, Lars
    J Jensen, Christian von Mering, Peer Bork Nucleic Acids Res. 2019 Jan 8;
    47(Database issue): D309–D314. doi: 10.1093/nar/gky1085 
```

Please, cite also the underlying algorithm used for the search step of eggNOG-mapper, and Prodigal if it was used for gene prediction:
```
[HMMER] Accelerated Profile HMM Searches. 
        Eddy SR. 2011. PLoS Comput. Biol. 7:e1002195.

[DIAMOND] Sensitive protein alignments at tree-of-life scale using DIAMOND.
          Buchfink B, Reuter K, Drost HG. 2021.
          Nature Methods 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x

[MMSEQS2] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
          Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988

[PRODIGAL] Prodigal: prokaryotic gene recognition and translation initiation site identification.
           Hyatt et al. 2010. BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119.

```