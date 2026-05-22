# Protein Feature Enrichment Score (PFES)

[![PyPI version](https://badge.fury.io/py/missense-pfes.svg)](https://pypi.org/project/missense-pfes/)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/broadinstitute/missense-pfes/blob/main/notebooks/G2P_PFES.ipynb)

Protein Feature Enrichment Score (PFES) is an interpretable framework for missense variant interpretation. Rather than outputting a single pathogenicity probability, PFES quantifies the degree to which a variant's protein characteristics statistically resemble known pathogenic or benign variation, and decomposes that signal into six feature attribute categories: physicochemical properties, 3D structure, domain/region annotations, functional sites, post-translational modifications, and protein-protein interactions.

This repository contains the data, analysis notebooks, and source code accompanying the manuscript: [*BioRxiv*](url)

---

## Repository Structure

```
missense-pfes/
├── data/                        # Pre-preprocessed protein feature annotations for case/control datasets
├── notebooks/                   # Analysis notebooks reproducing key figures and results, and live PFES computation (Colab)
├── results/                     # Output files - precomputed enrichment odds ratios and p-values, scored case and control datasets from batch pfes scorer, example output folder from PFES Colab
├── src/pfes/                    # PFES batch scorer (installable as a CLI tool)
└── glossary_protein_feature.md  # Definitions of all 103 protein features used for PFES computation
```

---
## Quick Start: Query a variant(s)

The easiest way to use PFES is through our **Google Colab notebook**, which requires no local installation: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/broadinstitute/missense-pfes/blob/main/notebooks/G2P_PFES.ipynb)

The notebook returns two outputs for any queried variant(s):


1. **PFES report** — overall PFES, PFES partitioning category with statistical significance, and a breakdown of contributions across six attribute categories with plain-language interpretation of each enriched feature

1. **Protein-wide mutational landscape** — a heatmap of PFES across all possible amino acid substitutions in the queried protein, decomposed by attribute category, showing where the variant of interest sits within the full substitution space


## PFES Batch Scorer

`src/pfes/scorer.py` computes PFES scores for a batch of missense variants from a TSV/CSV input file. It fetches protein feature data via the `g2papi` package and the enrichment table from this repository.

### Installation
```
pip install missense-pfes
```

Or install from source: 
```bash
git clone https://github.com/broadinstitute/missense-pfes.git
cd missense-pfes
pip install .
```

### Input format

The input file should be a TSV or CSV with the following columns:

| Column | Description |
|--------|-------------|
| `Gene` | HGNC gene symbol |
| `UniProt` | UniProt accession |
| `ResID` | Residue position |
| `RefAA` | Reference amino acid (single-letter) |
| `AltAA` | Alternate amino acid (single-letter) |

### Usage

```bash
usage: pfes [-h] -i INPUT -o OUTPUT [--rerun] [--log LOG] [--workers WORKERS]

Compute PFES scores for a batch of missense variants.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input TSV or CSV (columns: Gene, UniProt, ResID, RefAA, AltAA).
                        When --rerun is set, use the previous output file.
  -o OUTPUT, --output OUTPUT
                        Output TSV or CSV path (format inferred from extension).
  --rerun               Only process rows where PFES is NaN (fill missing scores).
  --log LOG             Path for error log (default: pfes_errors.log).
  --workers WORKERS     Number of parallel workers (default: 1).
```


```bash
# Basic usage
pfes -i variants.tsv -o scored_variants.tsv

# Parallel processing
pfes -i variants.tsv -o scored_variants.tsv --workers 4

# Specify a custom error log path
pfes -i variants.tsv -o scored_variants.tsv --log my_errors.log... 
```

Output retains the five input columns (`Gene`, `UniProt`, `ResID`, `RefAA`, `AltAA`) and appends seven score columns: `PFES`, `PFES_Physicochemical`, `PFES_Structure`, `PFES_Domain`, `PFES_Function`, `PFES_Modification`, and `PFES_PPI`.


### Handling failures and missing scores

Scoring requires fetching protein feature data from `g2papi` at runtime. Individual proteins or variants may fail due to missing data or network issues, in which case their PFES columns are left as `NaN` in the output. Details of each failure — including the gene, UniProt accession, residue, and error message — are written to a log file (`pfes_errors.log` by default).

To retry only the failed variants without re-running the entire dataset, pass the previous output file as input with `--rerun`:

```bash
# Fill in missing scores from a previous run
pfes -i scored_variants.tsv -o scored_variants.tsv --rerun
```

This fills in any rows where `PFES` is `NaN` and leaves successfully scored rows untouched.

## Data
Preprocessed protein feature annotations for the pathogenic (case) and control (benign + common population) variant datasets used in the paper. See `glossary_protein_feature.md` for definitions and sources of all 103 protein features.


| File | Description |
|------|-------------|
| `data/case_annotation.tsv` | Protein feature annotations for 39,695 pathogenic variants |
| `data/control_annotation.tsv` | Protein feature annotations for 130,832 control variants |


## Notebooks

| Notebook | Description |
|----------|-------------|
| `notebooks/compute_enrichment.ipynb` | Runs Fisher's exact test enrichment analysis within each PANTHER protein functional class and visualizes odds ratios. Generates `results/enrichment_OR_by_protein_class.csv`. |
| `notebooks/compute_partitioning_threshold.ipynb` | Computes empirical PFES distributions and derives partitioning thresholds (PF-Enriched/Neutral/Depleted) from the scored case and control datasets. |
| `notebooks/G2P_PFES.ipynb` |  Colab notebook for querying PFES scores and reports for individual variants of a gene of interest. <a href="https://colab.research.google.com/github/broadinstitute/missense-pfes/blob/main/notebooks/G2P_PFES_colab.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a> |

## Results

| File | Description |
|------|-------------|
| `results/case_pfes.tsv` | The scored (overall and per attribute) pathogenic variants, generated by the PFES batch scorer |
| `results/control_pfes.tsv` | The scored (overall and per attribute) control variants, generated by the PFES batch scorer |
| `results/pfes_pvalue_lookup.tsv` | Empirical p-value lookup table for variant partitioning, generated by `PFES_empirical_stats.ipynb`|
| `results/enrichment_OR_by_protein_class.csv` | Odds ratios from Fisher's exact tests for all 103 protein features across 20 PANTHER protein functional classes, generated by `enrichment_by_protein_class.ipynb` |


## Citation

If you use PFES in your work, please cite:

> Human Proteome-wide Mechanistic Interpretation of Missense Mutations through Protein Feature Enrichment Score *[BioRxiv DOI later]*

## Contact
Seulki Kwon — skwon@broadinstitute.org  
Broad Institute of MIT and Harvard
