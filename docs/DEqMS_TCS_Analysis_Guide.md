# DEqMS Analysis of Francisella TCS Proteins

## Overview

This document describes the complete workflow for differential expression analysis of Two-Component System (TCS) proteins in *Francisella* using DEqMS (Differential Expression analysis for proteomics with missing values).

### Key Features

- **PSM-dependent variance estimation** for improved statistical power
- **Multi-stress condition comparison** (oxidative, pH, heat, cold, iron limitation, etc.)
- **TCS protein-specific analysis** focusing on QseC/PmrA, KdpD/KdpE, and BfpR/BfpK systems
- **Publication-quality visualizations** including heatmaps, volcano plots, and correlation analyses

---

## Table of Contents

1. [Installation](#installation)
2. [Data Preparation](#data-preparation)
3. [Running the Analysis](#running-the-analysis)
4. [Output Files](#output-files)
5. [Biological Interpretation](#biological-interpretation)
6. [Troubleshooting](#troubleshooting)

---

## Installation

### Prerequisites

- R version ≥ 4.0.0
- RStudio (optional, but recommended)
- MaxQuant output files (proteinGroups.txt, evidence.txt)

### Install Required R Packages

```bash
# Run the installation script
Rscript scripts/install_r_packages.R
```

Or install manually in R:

```r
# CRAN packages
install.packages(c("dplyr", "tidyr", "ggplot2", "pheatmap",
                   "RColorBrewer", "gridExtra", "ggrepel"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "DEqMS", "ComplexHeatmap", "circlize"))
```

### Verify Installation

```r
# Load critical packages to verify
library(DEqMS)
library(limma)
library(ggplot2)
library(dplyr)
```

---

## Data Preparation

### Directory Structure

```
Francisella-ella-ella/
├── data/
│   ├── raw/
│   │   └── maxquant/
│   │       ├── proteinGroups.txt
│   │       └── evidence.txt
│   └── processed/
│       ├── PXD035145_proteinGroups.csv
│       └── PXD035145_PSM_counts.csv
├── src/R/deqms/
│   ├── deqms_tcs_analysis.R
│   ├── prepare_maxquant_data.R
│   └── visualize_tcs_results.R
├── results/
│   ├── tables/deqms/
│   └── figures/deqms/
└── scripts/
    └── install_r_packages.R
```

### Step 1: Prepare MaxQuant Data

Place your MaxQuant output files in `data/raw/maxquant/`:

- `proteinGroups.txt` - Main protein quantification file
- `evidence.txt` - Peptide-level data for PSM counting

Then run the preprocessing script:

```bash
Rscript src/R/deqms/prepare_maxquant_data.R
```

This script will:
- Filter contaminants and reverse hits
- Extract LFQ or intensity values
- Calculate PSM counts per protein
- Generate `PXD035145_proteinGroups.csv` and `PXD035145_PSM_counts.csv`

### Expected Data Format

**Protein Quantification Matrix:**
- Rows: Protein IDs (preferably gene names)
- Columns: Sample names
- Values: Raw intensities

**PSM Count Matrix:**
- Rows: Protein IDs (matching quantification matrix)
- Values: Number of peptide-spectrum matches

### Experimental Design

Update the `conditions` variable in `deqms_tcs_analysis.R` to match your experimental design:

```r
conditions <- factor(c(
  rep("Control", 3),
  rep("Oxidative", 3),
  rep("pH_stress", 3),
  rep("Heat_42C", 3),
  rep("Cold_25C", 3),
  rep("Iron_limit", 3),
  rep("Nutrient_starv", 3),
  rep("Osmotic", 3)
))
```

---

## Running the Analysis

### Complete Workflow

```bash
# 1. Install packages (one-time setup)
Rscript scripts/install_r_packages.R

# 2. Prepare MaxQuant data
Rscript src/R/deqms/prepare_maxquant_data.R

# 3. Run DEqMS analysis
Rscript src/R/deqms/deqms_tcs_analysis.R

# 4. Generate publication figures
Rscript src/R/deqms/visualize_tcs_results.R
```

### Step-by-Step Execution

#### Step 1: Main DEqMS Analysis

```bash
Rscript src/R/deqms/deqms_tcs_analysis.R
```

This performs:
1. Data import and preprocessing
2. Log2 transformation
3. Missing value filtering (≥20% valid values)
4. Median normalization
5. Linear modeling (limma)
6. DEqMS variance correction
7. Statistical testing for all contrasts
8. TCS protein extraction
9. Quality control plots

**Expected runtime:** 2-5 minutes for ~2000 proteins

#### Step 2: Visualization

```bash
Rscript src/R/deqms/visualize_tcs_results.R
```

Generates:
- Heatmaps of TCS expression
- Volcano plots for each contrast
- Expression profiles per TCS protein
- Correlation analyses
- Summary statistics

**Expected runtime:** 3-7 minutes

---

## Output Files

### Results Tables

Location: `results/tables/deqms/`

| File | Description |
|------|-------------|
| `deqms_all_results.csv` | Complete results for all proteins and contrasts |
| `TCS_deqms_all_results.csv` | All TCS protein results |
| `TCS_deqms_significant.csv` | Significant TCS changes only (|log2FC| > 1, p < 0.05) |
| `deqms_summary_statistics.csv` | Per-contrast summary statistics |

### Figures

Location: `results/figures/deqms/`

#### Quality Control Plots (`qc/`)

- `sample_distributions.pdf` - Boxplots of normalized intensities
- `missing_values.pdf` - Heatmap of missing value patterns
- `variance_psm_relationship.pdf` - DEqMS diagnostic plot

#### Main Figures

- `volcano_[contrast].pdf` - Volcano plot for each contrast highlighting TCS proteins
- `TCS_expression_heatmap.pdf` - Clustered heatmap of TCS fold changes

#### Publication Figures (`publication/`)

- `Fig1_TCS_heatmap_logFC.pdf` - Publication-ready heatmap
- `Fig2_volcano_panels.pdf` - Multi-panel volcano plots
- `Fig3_TCS_expression_profiles.pdf` - Individual TCS response profiles
- `Fig4_significance_summary.pdf` - Summary barplot
- `Fig5_pvalue_distribution.pdf` - P-value distributions
- `Fig6_TCS_correlation.pdf` - Correlation between TCS proteins
- `TableS1_TCS_DEqMS_results.csv` - Formatted supplementary table

---

## Biological Interpretation

### TCS Protein Systems in Francisella

#### QseC/PmrA System (Virulence Regulation)

**Function:**
- Master regulator of Francisella Pathogenicity Island (FPI)
- Controls expression of IglA, IglB, IglC, IglD, PdpA, etc.
- Essential for intracellular survival and replication

**Expected Responses:**
- ✓ **Upregulated** in oxidative stress (macrophage phagosome)
- ✓ **Upregulated** in iron limitation (host environment)
- ✓ **Upregulated** in nutrient starvation (ppGpp-mediated)
- ✗ No response to osmotic or temperature stress

**Locus Tags:**
- QseC: FTN_1617 (sensor kinase)
- PmrA: FTN_1465 (response regulator)

#### KdpD/KdpE System (Potassium Homeostasis)

**Function:**
- Senses potassium concentration and osmotic stress
- Regulates KdpABC high-affinity K+ transporter
- **NOTE:** Pseudogene in F. tularensis SchuS4 (functional in F. novicida)

**Expected Responses:**
- ✓ **Upregulated** in osmotic stress
- ✓ **Upregulated** in low K+ conditions
- ~ May respond to cold stress (membrane fluidity changes)

**Locus Tags:**
- KdpD: FTN_1714 (sensor kinase)
- KdpE: FTN_1715 (response regulator)

#### BfpR/BfpK System (Biofilm Formation)

**Function:**
- Regulates biofilm formation
- Inverse correlation with biofilm production
- Magnesium-responsive

**Expected Responses:**
- ✓ **Upregulated** under Mg²⁺ limitation
- ✓ **Upregulated** in response to antimicrobial peptides
- ~ Possible osmotic stress response

**Locus Tags:**
- BfpK: FTN_1453 (sensor kinase)
- BfpR: FTN_1452 (response regulator)

### Interpreting Results

#### Significant Upregulation (log2FC > 1, p < 0.05)

Indicates:
- TCS activation in response to specific stress
- Potential downstream regulatory cascade
- Biological relevance if matches expected pattern

#### Significant Downregulation (log2FC < -1, p < 0.05)

May indicate:
- Compensatory regulation
- Energy conservation under stress
- Indirect regulatory effects

#### No Significant Change

Could mean:
- TCS not responsive to this specific stress
- Insufficient statistical power (low PSM count)
- Post-translational regulation (not detected by proteomics)

### Cross-Referencing with FPI Proteins

Check if TCS changes correlate with FPI protein expression:

```r
# Example: Check if IglC is upregulated with QseC/PmrA
fpi_genes <- c("IglA", "IglB", "IglC", "IglD", "PdpA", "PdpB")
fpi_results <- all_results %>%
  filter(gene %in% fpi_genes)
```

---

## Statistical Considerations

### DEqMS vs. Standard limma

**DEqMS advantages:**
- Accounts for heterogeneous quantification quality
- More accurate variance estimation
- Improved statistical power for low-PSM proteins

**When to use DEqMS:**
- Proteomics data with missing values
- Variable PSM counts across proteins
- Need for maximum statistical power

### Significance Thresholds

**Default thresholds:**
- Log2 fold change: |log2FC| > 1 (2-fold change)
- Adjusted p-value: < 0.05 (FDR-corrected)

**Adjusting thresholds:**

```r
# More stringent
LOG2FC_THRESHOLD <- 1.5  # 3-fold change
PVAL_THRESHOLD <- 0.01   # 1% FDR

# More permissive
LOG2FC_THRESHOLD <- 0.58  # 1.5-fold change
PVAL_THRESHOLD <- 0.1     # 10% FDR
```

### Multiple Testing Correction

DEqMS uses Benjamini-Hochberg FDR correction by default:

- `sca.P.Value` - Nominal p-value
- `sca.adj.pval` - FDR-adjusted p-value (**use this for significance**)

---

## Troubleshooting

### Common Issues

#### 1. "No LFQ intensity columns found"

**Cause:** MaxQuant was not run with LFQ option

**Solution:**
```r
# In prepare_maxquant_data.R, it will automatically fall back to:
extract_intensities(protein_groups, intensity_type = "Intensity")
```

#### 2. "PSM counts not matching protein IDs"

**Cause:** Gene name / Protein ID mismatch between proteinGroups and evidence

**Solution:** Script automatically uses peptide counts as proxy if evidence.txt is missing

#### 3. "Too many missing values"

**Cause:** Stringent filtering removes most proteins

**Solution:** Adjust `MIN_VALID_RATIO`:
```r
MIN_VALID_RATIO <- 0.1  # Allow 10% valid values (less stringent)
```

#### 4. "No TCS proteins found in data"

**Cause:** Locus tags don't match your organism

**Solution:** Update `TCS_GENES` mapping:
```r
TCS_GENES <- c(
  "FTT_1557c" = "QseC",  # Example for F. tularensis SchuS4
  # ... add your organism's locus tags
)
```

#### 5. "Package installation failed"

**Cause:** Missing system dependencies

**Ubuntu/Debian solution:**
```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

**macOS solution:**
```bash
brew install openssl libxml2
```

### Getting Help

For issues specific to:
- **DEqMS methodology:** [DEqMS paper (MCP 2020)](https://doi.org/10.1074/mcp.TIR119.001646)
- **limma usage:** [limma user guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
- **Francisella TCS:** Check references in main analysis documentation

---

## Advanced Customization

### Adding Custom Contrasts

```r
# Example: Compare stress conditions to each other
contrasts_matrix <- makeContrasts(
  Heat_vs_Cold = Heat_42C - Cold_25C,
  Oxidative_vs_pH = Oxidative - pH_stress,
  # ... additional contrasts
  levels = design
)
```

### Analyzing Additional Protein Sets

```r
# Example: Analyze all response regulators
RR_GENES <- c("FTN_1465", "FTN_1715", "FTN_1452", "FTN_0527", ...)

rr_results <- all_results %>%
  filter(gene %in% RR_GENES)
```

### Exporting for Pathway Analysis

```r
# Export significant proteins for DAVID/STRING analysis
sig_proteins <- all_results %>%
  filter(
    contrast == "Oxidative_vs_Control",
    abs(logFC) > 1,
    sca.adj.pval < 0.05
  ) %>%
  pull(gene)

write.table(sig_proteins, "sig_proteins_oxidative.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
```

---

## Citation

If you use this workflow, please cite:

**DEqMS:**
> Zhu Y, Orre LM, Tran YZ, et al. DEqMS: A Method for Accurate Variance Estimation in Differential Protein Expression Analysis. *Mol Cell Proteomics*. 2020;19(6):1047-1057. doi:10.1074/mcp.TIR119.001646

**limma:**
> Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47. doi:10.1093/nar/gkv007

---

## License

This workflow is provided as-is for academic research purposes.

---

## Contact

For questions or issues with this workflow, please open an issue on the project repository.

Last updated: 2025
