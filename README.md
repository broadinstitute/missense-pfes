# README

## Protein Feature Enrichment Score (PFES) Script

This Python script generates heatmaps and analyses mutation scores for protein features based on enrichment scores provided in data files. The heatmaps and mutation analyses are saved as output files for further interpretation.



# 1.  `get_pfes_score.py`

## Requirements

### Libraries:

* `python 3`
* `pandas`
* `numpy`
* `matplotlib`
* `seaborn`
* `argparse`
* `requests`
* `scipy`

### External Tools:

* `gsutil` for accessing files stored in Google Cloud Storage.
  ```
  pip install google-cloud-storage
  ```

## Usage 

### Command Line Arguments

This script accepts the following arguments:

#### Required (Mutually Exclusive):

* `-gene`: The gene name to process.
* `-uniprot`: The UniProt ID to process.

#### Optional:

* `-mutation`: A single mutation (e.g., `M1A`) to analyze.
* `-mutations`: A file containing a list of mutations (one per line).

### Example Command:

```bash
python script.py -gene BID -mutation H42K
```

Or with a mutations file:

```bash
python script.py -uniprot P55957 -mutations mutations.txt
```

---

## Input

### Data Files:

1. **UniProt Metadata** :

* File: `uniprot_metadata_2024_04.tsv`
* Format: Tab-delimited with required columns: `UniprotKB_Entry`.

1. **Gene Metadata** :

* File: `gene_metadata_2024_04.tsv`
* Format: Tab-delimited with required columns: `HGNC_symbol` and `UniprotKB_Entry`.

1. **PLP Likelihood Data** :

* File: `plp_likelihood.tsv`
* Format: Tab-delimited with columns `Center` and `PLP dominance`.

Above files are provided in the folder `./files/`

### Mutation File (Optional):

* Text file containing one mutation per line (e.g., `M1A`).
* Example file is given as `mutations.txt`

## Output

### Generated Files:

1. **Heatmaps** :

* File: `heatmap_{key}_{gene_name}_{uid}.png`
* Location: `../data_{gene_name}_{uid}`
* Description: Heatmaps of enrichment scores for various features.

1. **Mutation Analysis CSV** :

* File: `mutations_PFES_{gene_name}_{uid}.csv`
* Location: `../data_{gene_name}_{uid}`
* Description: Contains detailed mutation scores, pathogenicity likelihood, and decomposed PFE scores.

### Example CSV Columns:

* `Mutation`: The input mutation.
* `Total PFE Score`: The overall score for the mutation.
* `Pathogenicity Likelihood`: Likelihood of pathogenicity (in %) derived from precalculated statistical analysis.
* `Physicochemical/Function/Domain/Modification/Structure/PPI PFE Score`: Decomposed score for each mechanisms features.
* `List of features;Odd ratio`: Features and their respective odd ratios.

# 2. get_pfes_score_interactive_notebook.ipynb

This notebook provides interactive exploration of the heatmap and relavant features for a given gene and mutations. The output is largely same as above. It gives additional information on

* Corresponding protein class and its enriched feature list
* Pathogenicity likelihood curve for a given mutation
* A list of features for a given feature category (physicochemical, function, domain, modifications, structure, and PPI)
