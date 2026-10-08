# KlaAtlas

KlaAtlas is a cross-species resource for exploring lysine lactylation at protein coordinates.

**Visit [www.kla-atlas.com](https://www.kla-atlas.com/)** to search proteins and sites, browse annotation tracks, trace study evidence, build candidate shortlists and export selected records.

## Database and sources

`database/klaatlas.sqlite` contains 142,204 Kla sites, 35,613 protein–taxon entries and 334,986 evidence records across 23 taxa and 115 publications. It also contains the reference sequences, analysis inputs and statistical results.

`database/DOI_sources.tsv` lists the source publications and DOI links.

```sh
git lfs install
git clone https://github.com/MforMegaptera/KlaAtlas.git
cd KlaAtlas
git lfs pull
```

## Analyses and figures

`figures/analysis/` contains the analysis code, `figures/code/` contains the plotting code, and `figures/figures/` contains the figures and visual source files.

Use Python 3.11, a C99 compiler (`cc`), and R 4.4.3 with survival 3.8-3.

```sh
python -m pip install -r figures/requirements.txt

# Analysis and plotting
python figures/run.py --stage all --figure all \
  --database database/klaatlas.sqlite --output results

# Analysis for one figure
python figures/run.py --stage analysis --figure 4 \
  --database database/klaatlas.sqlite --output results

# Plot the analysis results
python figures/run.py --stage plot --figure 4 \
  --database database/klaatlas.sqlite --output results --analysis-dir results/analysis

# Plot from the database
python figures/run.py --stage plot --figure all \
  --database database/klaatlas.sqlite --output results
```

Software: [MIT](LICENSE). Reference sequences: The UniProt Consortium ([CC BY 4.0](https://www.uniprot.org/help/license)). KEGG annotations: [source terms](https://www.kegg.jp/kegg/legal.html).
