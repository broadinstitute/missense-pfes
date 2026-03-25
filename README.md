# Protein Feature Enrichment Score (PFES)

> **Human Proteome-wide Mechanistic Interpretation of Missense Mutations through Structural and Functional Characterization** 
>
> Seulki Kwon, Jordan Safer, Sumaiya Iqbal (Broad Institute of MIT and Harvard) 

---

## Overview

Accurately interpreting missense variants remains a central challenge in clinical genetics. Most computational variant effect predictors achieve strong classification performance but function as black boxes, providing pathogenicity probabilities without revealing the protein characteristics that drive them.

**P**rotein **F**eature **E**nrichment **S**core (**PFES**) is an variant interpretation framework that asks a different question: *which protein characteristics are implicated in a variant's molecular profile, and how does that profile relate to known pathogenic signatures?*

Rather than optimizing for classification accuracy, PFES quantifies the degree to which a missense variant's protein context resembles known pathogenic or benign variants, and decomposes that signal into six interpretable feature attribute categories. This makes PFES a complementary resource to existing variant effect predictors, adding a transparent mechanistic layer to variant interpretation workflows.

---

## The PFES Framework

### Core principle

The framework is grounded in the idea that **molecular context determines consequence**: pathogenicity emerges from the intersection of *what* is being changed (physicochemical properties), *where* the change occurs (structural and functional context), and *why* that location matters (protein class-specific constraints).

### Score computation

PFES is computed as the sum of log odds ratios across protein features showing statistically significant enrichment in pathogenic versus benign/population variants:

$$\text{PFES} = \sum_{i \in \text{significant}} \log(\text{OR}_i)$$

Features enriched in pathogenic variants contribute positive values; features enriched in benign/population variants contribute negative values. Statistical enrichment is assessed using Fisher's exact tests with Benjamini-Hochberg FDR correction (threshold: p < 0.01), performed separately within each of 19 PANTHER protein functional classes to capture class-specific biological constraints.

PFES spans **103 protein features** across six attribute categories:

| Attribute | Examples |
|---|---|
| Physicochemical | Grantham distance, amino acid group |
| Structure | Relative surface area, secondary structure, intramolecular interactions, AlphaFold pLDDT |
| Domain/Region | Transmembrane regions, disordered regions, signal peptides |
| Function | Active sites, binding sites, DNA-binding regions |
| Modification | Phosphorylation, acetylation, methylation, lipidation sites |
| PPI | Intermolecular hydrogen bonds, van der Waals contacts, disulfide bonds |

### Variant partitioning

Variants are partitioned into three categories based on statistical deviation from empirical PFES distributions of clinically annotated pathogenic and benign variants:

- **PF-Enriched**: PFES significantly deviates from the benign/population distribution (p < 0.05), indicating a protein feature profile statistically consistent with pathogenic variation
- **PF-Depleted**: PFES significantly deviates from the pathogenic distribution (p < 0.05), indicating a protein feature profile statistically consistent with benign variation
- **PF-Neutral**: PFES is statistically consistent with both distributions, lacking a distinctive protein feature signal in either direction

This is **mechanistic partitioning**, not binary pathogenicity classification. PF-neutral status should not be interpreted as evidence of benignity; it reflects that the variant lacks a distinctive protein feature profile detectable by this framework, and other lines of evidence remain essential.

---

## Repository Contents

```
pfes/
├── README.md
├── enrichment/               # Scripts for PFES calculation
│   ├── ...
├── notebooks/
│   └── pfes_query.ipynb      # Google Colab notebook for variant query
├── data/
│   └── ...                   # Precomputed PFES for all human proteome missense variants
└── ...
```

---

## Quick Start: Query a Variant

The easiest way to use PFES is through our **Google Colab notebook**, which requires no local installation:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](TBU)

The notebook returns two outputs for any queried variant:

1. **Variant summary report** — overall PFES, PF category with statistical significance, and a breakdown of contributions across six attribute categories with plain-language interpretation of each enriched feature
2. **Protein-wide mutational landscape** — a heatmap of PFES across all possible amino acid substitutions in the queried protein, decomposed by attribute category, showing where the variant of interest sits within the full substitution space

---

## Interpreting PFES Output

PFES is most informative when used alongside existing variant effect predictors rather than in isolation.

- A variant that is **PF-Enriched** and receives high pathogenicity scores from tools like AlphaMissense or REVEL carries a richer, more interpretable evidence profile: the existing tools address whether the variant is likely pathogenic, while PFES addresses which specific protein characteristics are implicated.
- A variant scored as likely pathogenic by other tools but **PF-Neutral** may involve a mechanism outside the scope of currently annotated protein features (e.g., gain-of-function through novel interactions, subtle allosteric effects), warranting closer investigation.
- **PF-Depleted** status can support benign reclassification as one molecular line of evidence among others, particularly when corroborated by population frequency and clinical data.

Because PFES is derived from protein feature annotations rather than trained on labeled variant data, it provides evidence that is genuinely orthogonal to frequency-based, segregation-based, and co-occurrence-based evidence lines used in ACMG/AMP classification workflows.

---

## Limitations

- PFES is entirely dependent on existing protein annotations (PDB, AlphaFold, UniProtKB, PhosphoSitePlus, PANTHER). Poorly characterized proteins will yield lower-resolution profiles.
- The framework detects disruption of protein features statistically enriched at the proteome or protein class level. Variants acting through highly protein-specific mechanisms, subtle allosteric perturbations, or gain-of-function effects involving unannotated interactions will often be PF-Neutral.
- Because PFES is primarily tied to protein features at the reference amino acid position, it is better suited to characterizing constraint landscapes across positions within a protein than distinguishing between individual substitutions at a single site.
- Enrichment patterns are derived from ClinVar and gnomAD, which are not uniformly distributed across the proteome. Protein classes underrepresented in clinical databases may have less robust enrichment estimates.

---

## Citation

If you use PFES in your work, please cite:

> Kwon S, Safer J, Iqbal S. *Interpretable Protein Feature Enrichment Analysis for Human Proteome-wide Missense Variant Classification.* (2025) [journal TBU]

---

## Contact

For questions or feedback, please open a GitHub issue or contact [sumaiya@broadinstitute.org](mailto:sumaiya@broadinstitute.org).
