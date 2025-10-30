# Francisella Two-Component Systems Proteomics Analysis

## Overview
Computational framework for analyzing two-component systems in Francisella species using proteomics data.

## Project Status
✅ Working - Demo version with sample data

## Author
- Brittany Moretz (bmoretz@gmu.edu.edu)

## Description
This repository contains code and analysis for identifying TCS-regulated secretion networks in Francisella novicida. The analysis pipeline processes proteomics data to:

- Identify two-component system proteins (QseB, PmrA, BfpR, KdpE, QseC, KdpD, etc.)
- Analyze phosphorylation sites on TCS components
- Predict secreted proteins using sequence-based methods
- Perform differential expression analysis
- Integrate results to find TCS-regulated secreted proteins

## Quick Start

### Installation

1. Clone this repository:
```bash
git clone https://github.com/bam-builds/francisella-ella-ella.git
cd francisella-ella-ella
```

2. Install required Python packages:
```bash
pip install -r requirements.txt
```

3. (Optional) Install R packages for advanced visualizations:
```bash
Rscript src/R/install_packages.R
```

### Running the Analysis

#### Python Analysis Pipeline

Run the main analysis script:
```bash
python3 analysis.py
```

This will:
- Generate sample proteomics data for demonstration
- Analyze TCS proteins and phosphorylation sites
- Predict secreted proteins
- Perform differential expression analysis
- Create visualizations in `results/figures/`
- Export results to `results/tables/`

#### R Visualization Pipeline

Generate publication-quality figures:
```bash
Rscript src/R/tcs_proteomics_visualization.R
```

This will create:
- 7 comprehensive publication-quality figures
- 1 composite multi-panel figure
- Statistical summary tables
- Figure legends for manuscripts

See `src/R/README.md` for detailed documentation.

### Output Files

After running the analysis, you'll find:

**Python Analysis Figures:**
- `results/figures/phosphosite_analysis.png` - Phosphorylation site analysis
- `results/figures/volcano_plot.png` - Differential expression volcano plot

**R Visualization Figures:**
- `figures/Fig1_TCS_Regulation_Heatmap.pdf` - Comprehensive TCS heatmap
- `figures/Fig2_Mutant_Comparison.pdf` - WT vs mutant strain comparison
- `figures/Fig3_Temporal_Dynamics.pdf` - Time-course infection analysis
- `figures/Fig4_Volcano_Plots.pdf` - Multi-condition volcano plots
- `figures/Fig5_Regulatory_Network.pdf` - TCS regulatory network
- `figures/Fig6_Species_Comparison.pdf` - Cross-species analysis
- `figures/Fig7_Summary_Statistics.pdf` - Statistical summaries
- `figures/Fig_Composite_TCS_Analysis.pdf` - Publication composite figure

**Data Tables:**
- `results/tables/tcs_phosphorylation_sites.csv` - TCS phosphorylation sites
- `results/tables/tcs_regulated_secreted_proteins.csv` - Integrated results
- `results/tables/differential_expression.csv` - DE analysis results
- `figures/Table_S1_TCS_Significant_Changes.csv` - Supplementary TCS changes

## Project Structure

```
francisella-ella-ella/
├── analysis.py                 # Main Python analysis script
├── requirements.txt            # Python dependencies
├── src/
│   ├── python/
│   │   ├── data_retrieval.py   # ProteomeXchange data download
│   │   ├── secretion_prediction.py  # Secretion signal prediction
│   │   └── ptm_analysis.py     # PTM analysis tools
│   └── R/
│       ├── tcs_proteomics_visualization.R  # Main R visualization script
│       ├── install_packages.R  # R package installer
│       └── README.md           # R-specific documentation
├── data/
│   ├── raw/                    # Raw data files
│   ├── processed/              # Processed data
│   └── metadata/               # Metadata files
├── results/
│   ├── figures/                # Python-generated plots
│   └── tables/                 # Analysis results
└── figures/                    # R-generated publication figures
```

## Key Features

### Python Analysis Pipeline
- **TCS Protein Identification**: Identifies response regulators, sensor kinases, and FPI genes
- **Phosphoproteomics**: Analyzes phosphorylation sites with localization probability filtering
- **Secretion Prediction**: Predicts proteins with secretion signals
- **Differential Expression**: Identifies significantly changed proteins
- **Integration**: Finds TCS-regulated secreted proteins

### R Visualization Pipeline
- **Publication-Quality Figures**: 7+ comprehensive figures ready for publication
- **Complex Heatmaps**: Hierarchical clustering with statistical annotations
- **Network Visualization**: TCS regulatory cascade networks
- **Temporal Analysis**: Time-course infection dynamics
- **Cross-Species Comparison**: Functional capacity across Francisella strains
- **Statistical Summaries**: Comprehensive expression pattern analysis

## Requirements

### Python (Analysis)
- Python 3.7+
- pandas
- numpy
- matplotlib
- seaborn
- requests

See `requirements.txt` for full Python dependencies.

### R (Visualization)
- R 4.0+
- tidyverse
- ComplexHeatmap
- ggplot2 and extensions (ggrepel, cowplot, etc.)
- Network visualization packages

Run `Rscript src/R/install_packages.R` to install all R dependencies.

## Future Development

- Integration with real ProteomeXchange datasets (PXD009225, PXD035145, etc.)
- Advanced machine learning models for secretion prediction
- Protein interaction network analysis
- Time-series analysis for infection studies

