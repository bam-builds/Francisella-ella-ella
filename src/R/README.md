# TCS Proteomics Visualization (R)

## Overview

This directory contains R scripts for creating publication-quality visualizations of Two-Component System (TCS) proteomics data from Francisella species.

## Features

The visualization pipeline generates 7 comprehensive figures plus a composite publication figure:

1. **Figure 1**: Comprehensive TCS Regulation Heatmap
   - Hierarchical clustering of TCS proteins across stress conditions
   - Statistical significance annotations
   - Condition and protein type annotations

2. **Figure 2**: Wild-Type vs Mutant Strain Comparison
   - Comparative barplots of TCS expression
   - Multiple isogenic deletion strains (ΔmglA, ΔqseC, ΔqseB)
   - Grouped by protein category (TCS, FPI, Stress Response)

3. **Figure 3**: Temporal Dynamics During Infection
   - Time-course analysis (0-24 hours post-infection)
   - TCS, FPI, and global regulator expression patterns
   - Annotated infection phases (phagosomal vs cytosolic)

4. **Figure 4**: Volcano Plot Grid
   - Four key conditions: Oxidative stress, pH stress, Iron limitation, Macrophage infection
   - Differential expression with statistical thresholds
   - Labeled key proteins

5. **Figure 5**: Regulatory Cascade Network
   - Network visualization of TCS regulatory pathways
   - Node coloring by expression level
   - Edge weight by interaction strength

6. **Figure 6**: Species Comparison
   - Cross-species TCS functional capacity
   - F. novicida vs pathogenic strains (SchuS4, LVS)
   - Pseudogene/gene loss annotation

7. **Figure 7**: Summary Statistics
   - Bar plots of regulated protein counts
   - Line plots of average expression changes
   - TCS vs FPI comparison

8. **Composite Figure**: Multi-panel publication figure combining key visualizations

## Installation

### Prerequisites

- R version 4.0 or higher
- Sufficient memory for complex heatmaps (recommended: 8GB+ RAM)

### Install Required Packages

Run the installation script to install all required R packages:

```bash
Rscript src/R/install_packages.R
```

This will install:

**CRAN Packages:**
- tidyverse, dplyr, ggplot2, tidyr
- viridis, RColorBrewer, ggrepel
- cowplot, ggsci, scales, gridExtra
- network, ggnetwork, pheatmap, circlize

**Bioconductor Packages:**
- ComplexHeatmap

## Usage

### Basic Usage (Mock Data)

The script includes mock data generation for demonstration:

```bash
cd /path/to/Francisella-ella-ella
Rscript src/R/tcs_proteomics_visualization.R
```

This will create a `figures/` directory in your current working directory with all output PDFs.

### Using Real Data

To use your own proteomics data, modify the `generate_mock_tcs_data()` function or replace it with a data loading function:

```r
# Replace mock data generation with actual data loading
load_real_tcs_data <- function(file_path) {
  # Load your proteomics data
  expression_data <- read.csv(file_path)

  # Process into required format:
  # - expression: matrix with proteins as rows, conditions as columns
  # - p_values: matrix of p-values
  # - mutant_data: data frame for strain comparisons

  return(list(
    expression = expression_matrix,
    p_values = p_value_matrix,
    mutant_data = mutant_comparisons,
    proteins = protein_names
  ))
}

# In main execution section:
tcs_data <- load_real_tcs_data("data/processed/tcs_expression.csv")
```

### Expected Data Format

The visualization functions expect data in the following format:

**Expression Matrix:**
```
                Control  Oxidative_2h  Oxidative_8h  ...
QseB/PmrA       0.00     1.45          -0.89        ...
QseC            0.00     2.13           0.52        ...
...
```

**P-values Matrix:**
```
                Oxidative_2h  Oxidative_8h  ...
QseB/PmrA       0.001         0.045        ...
QseC            0.0001        0.023        ...
...
```

**Mutant Comparison Data Frame:**
```
protein    strain    expression  category
QseB/PmrA  WT        0.00        TCS
QseB/PmrA  ΔmglA    -2.50        TCS
...
```

## Output Files

All figures are saved as PDF files in the `figures/` directory:

- `Fig1_TCS_Regulation_Heatmap.pdf` (12" × 8")
- `Fig2_Mutant_Comparison.pdf` (14" × 6")
- `Fig3_Temporal_Dynamics.pdf` (10" × 10")
- `Fig4_Volcano_Plots.pdf` (12" × 10")
- `Fig5_Regulatory_Network.pdf` (10" × 8")
- `Fig6_Species_Comparison.pdf` (12" × 8")
- `Fig7_Summary_Statistics.pdf` (10" × 10")
- `Fig_Composite_TCS_Analysis.pdf` (16" × 14")
- `Table_S1_TCS_Significant_Changes.csv`
- `Figure_Legends.txt`

## Customization

### Color Schemes

The script uses colorblind-friendly palettes. To customize:

```r
# Heatmap colors (Figure 1)
col = colorRamp2(c(-4, -2, 0, 2, 4),
                 c("#2166AC", "#4393C3", "white", "#F4A582", "#B2182B"))

# Strain colors (Figure 2)
scale_fill_manual(values = c("WT" = "#4DAF4A", "ΔmglA" = "#E41A1C", ...))
```

### Figure Dimensions

Adjust PDF dimensions in the main execution section:

```r
pdf("figures/Fig1_TCS_Regulation_Heatmap.pdf", width = 12, height = 8)
```

### Statistical Thresholds

Modify significance thresholds in volcano plots:

```r
significance = case_when(
  padj < 0.05 & abs(log2FC) > 1.5 ~ "Significant",  # Change thresholds here
  TRUE ~ "Not Significant"
)
```

## Troubleshooting

### Memory Issues

For large datasets, increase R memory limit:

```r
# At the beginning of the script
memory.limit(size = 16000)  # 16GB on Windows
```

### Font Issues

If you encounter font warnings, install additional fonts or use base fonts:

```r
theme_set(theme_minimal(base_family = "sans"))
```

### Missing Packages

If automatic installation fails, install packages manually:

```r
install.packages("tidyverse")
BiocManager::install("ComplexHeatmap")
```

## Integration with Python Pipeline

This R visualization pipeline complements the Python analysis pipeline:

1. **Python analysis** (`analysis.py`):
   - Data retrieval from ProteomeXchange
   - Phosphoproteomics analysis
   - Differential expression analysis
   - Secretion prediction

2. **R visualization** (`tcs_proteomics_visualization.R`):
   - Publication-quality figures
   - Complex heatmaps
   - Network visualizations
   - Statistical summaries

### Workflow

```bash
# Step 1: Run Python analysis
python3 analysis.py

# Step 2: Generate R visualizations
Rscript src/R/tcs_proteomics_visualization.R
```

## Citation

If you use this visualization pipeline in your research, please cite:

```
Moretz, B. (2024). Francisella Two-Component Systems Proteomics Analysis.
GitHub repository: https://github.com/bam-builds/francisella-ella-ella
```

## License

This project is part of the Francisella TCS proteomics analysis framework.

## Contact

For questions or issues:
- Email: bmoretz@gmu.edu.edu
- GitHub: https://github.com/bam-builds/francisella-ella-ella/issues

## References

Key TCS proteins analyzed:
- **QseB/PmrA**: Response regulator (orphan system)
- **QseC**: Sensor histidine kinase
- **KdpD/KdpE**: Potassium-sensing two-component system
- **BfpR/BfpK**: Biofilm-related two-component system
- **MglA/SspA**: Global transcriptional regulators

FPI (Francisella Pathogenicity Island) genes:
- **IglA-D**: Intracellular growth locus proteins
- **PdpA/B**: Pathogenicity determinant proteins
- **VgrG**: Type VI secretion effector
