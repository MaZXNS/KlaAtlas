# KlaAtlas: Data and Analysis Figures

KlaAtlas is a cross-species resource for exploring lysine lactylation at protein coordinates. Visit **[www.kla-atlas.com](https://www.kla-atlas.com/)** to explore the atlas.

This package uses the November 4, 2026 database version, `ROOT-MIN10-CONSERVATIVE-GROUP-20261104`, with a literature search cutoff of October 24, 2026. It contains **212,444 sites, 568,654 direct evidence records, 46,280 protein–taxon entries and 247 source tables**. The SQLite database provides curated sites, evidence, sources, reference sequences and figure inputs derived from public result tables.

- [Database and query guide](database/README.md)
- [Analysis figures](figures/)
- [Data dictionary](database/DATA_DICTIONARY.tsv)
- [File checksums](SHA256SUMS)

## Getting started

Use Python 3.11 and install the Arial font for plotting.

```sh
git lfs install
git clone https://github.com/MaZXNS/KlaAtlas.git
cd KlaAtlas
git lfs pull
python verify_package.py
python -m pip install -r requirements.txt
python database/restore_database.py
python verify_package.py --database database/klaatlas.sqlite
python make_figures.py --database database/klaatlas.sqlite --output generated_figures
```

## Analysis figures

| Figure | Content |
| --- | --- |
| Figure 1b–f | Atlas coverage, source years, coordinate agreement and cumulative coverage. |
| Figure 2 | Site coverage, positional distribution, clustering, peptide and localization statistics. |
| Figure 3 | Sequence features, symmetry and histone analyses. |
| Figure 4 | Subcellular localization, structure and functional enrichment. |

Select a figure with `--figure 1`, `2`, `3` or `4`, or use `--figure all` to generate all analysis figures. The output `run_metadata.json` records input hashes, run parameters and software versions.

## Data

`sites` and `evidence` contain sites and their source observations. `proteins`, `sources`, `publications` and `sequences` provide protein, source and sequence information. `figure_data` stores the analysis figure inputs. Evidence can be traced through `source_file_sha256`, `source_table_id`, `source_sheet_or_table` and `source_row_number`.
