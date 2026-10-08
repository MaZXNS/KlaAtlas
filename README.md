# KlaAtlas

KlaAtlas is a cross-species resource for exploring lysine lactylation at protein coordinates.

**Visit [www.kla-atlas.com](https://www.kla-atlas.com/)** to search proteins and sites, browse sequence-aligned annotation tracks, trace study evidence, build candidate shortlists and export selected records.

## Database and sources

`database/klaatlas.sqlite` contains 142,204 unique Kla sites, 35,613 protein–taxon entries and 334,986 source-evidence records across 23 taxa and 115 publications. Exact reference sequences and source relationships are retained. `database/DOI_sources.tsv` links publication identifiers to the 115 source DOIs.

The same database stores the frozen data needed for the analyses: comparison populations, covariates, eligibility and censoring flags, reference-annotation memberships, plotting inputs and complete statistical test results. Compressed input relations are extracted automatically to temporary working files. No original literature PDFs, publisher supplementary files or raw mass-spectrometry files are distributed.

A site is identified by NCBI taxon ID, exact UniProt accession and 1-based lysine position on its reference sequence. Evidence-record counts are source-support counts, rather than numbers of independent experiments. Unreported lysines and missing annotations are not experimentally established negative results.

Download the database with Git LFS:

```sh
git lfs install
git clone https://github.com/MforMegaptera/KlaAtlas.git
cd KlaAtlas
git lfs pull
```

## Analyses and figures

`figures/` contains the statistical analysis and plotting code, with one entry point for analysis, drawing or both. Analyses use the frozen processed inputs and reference features in the database; they do not refresh online annotations or reprocess raw mass-spectrometry files.

Use Python 3.11, a C99 compiler (`cc`), and R 4.4.3 with survival 3.8-3 for the original Efron conditional logistic model.

```sh
python -m pip install --upgrade pip
python -m pip install -r figures/requirements.txt

# Run all statistical analyses and generate the figures.
python figures/run.py --stage all --figure all \
  --database database/klaatlas.sqlite --output results

# Run one figure's analyses, then draw from those new results.
python figures/run.py --stage analysis --figure 4 \
  --database database/klaatlas.sqlite --output results
python figures/run.py --stage plot --figure 4 \
  --database database/klaatlas.sqlite --output results --analysis-dir results/analysis

# Draw the frozen results without refitting the statistical models.
python figures/run.py --stage plot --figure all \
  --database database/klaatlas.sqlite --output frozen_figures
```

| Figure | Analysis inputs and retained methods |
|---|---|
| 1 | Atlas membership, study/source blocks and exact-reference K opportunities; resource counts, 200 randomized source orders and source-pair positional comparisons. Panel a is a fixed curation workflow. |
| 2 | Protein K opportunities and observed/theoretical peptide model inputs; original conditional models, 999 fixed-burden placements, 10,000 placements per protein and 1,000 protein bootstrap replicates. |
| 3 | Same-reference sequence windows, matched strata and histone/comparison populations; all 400 position tests, 200 paired contrasts and the original histone coefficient families. |
| 4 | Localization bindings, structural lysine features and 13 measured proteomes with abundance/K covariates; compartment models, pLDDT ≥90/70 Efron models and all 15,783 functional tests with BH correction separately within 52 study × namespace families. |
| 5 | Fixed interface and protein-track data; the website view is rendered from its recorded source material. |
| S1 | Fixed source diagram of study methods; supplied as a project figure with its editable visual source. |

Original filtering rules, missing-data handling, model formulas, correction families and random seeds are retained. Drawing frozen results is separate from refitting models. The atom-level construction of the supplied structural features is outside this processed-input workflow. Failed estimates and unavailable annotations remain explicit.

Software license: [MIT](LICENSE). UniProt-derived references retain attribution to The UniProt Consortium ([CC BY 4.0](https://www.uniprot.org/help/license)); other upstream study and annotation terms remain applicable. KEGG-derived memberships retain [KEGG's source terms](https://www.kegg.jp/kegg/legal.html), including its academic-use restrictions, and are not relicensed under MIT.
