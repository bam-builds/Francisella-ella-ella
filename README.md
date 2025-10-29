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
git clone https://github.com/bam-builds/francisella-tcs-analysis.git
cd francisella-tcs-analysis
```

2. Install required Python packages:
```bash
pip install -r requirements.txt
```

### Running the Analysis

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

### Output Files

After running the analysis, you'll find:

**Figures:**
- `results/figures/phosphosite_analysis.png` - Phosphorylation site analysis
- `results/figures/volcano_plot.png` - Differential expression volcano plot

**Data Tables:**
- `results/tables/tcs_phosphorylation_sites.csv` - TCS phosphorylation sites
- `results/tables/tcs_regulated_secreted_proteins.csv` - Integrated results
- `results/tables/differential_expression.csv` - DE analysis results

## Project Structure

```
francisella-tcs-analysis/
├── analysis.py                 # Main analysis script
├── requirements.txt            # Python dependencies
├── src/
│   └── python/
│       ├── data_retrieval.py   # ProteomeXchange data download
│       ├── secretion_prediction.py  # Secretion signal prediction
│       └── ptm_analysis.py     # PTM analysis tools
├── data/
│   ├── raw/                    # Raw data files
│   ├── processed/              # Processed data
│   └── metadata/               # Metadata files
└── results/
    ├── figures/                # Generated plots
    └── tables/                 # Analysis results
```

## Key Features

- **TCS Protein Identification**: Identifies response regulators, sensor kinases, and FPI genes
- **Phosphoproteomics**: Analyzes phosphorylation sites with localization probability filtering
- **Secretion Prediction**: Predicts proteins with secretion signals
- **Differential Expression**: Identifies significantly changed proteins
- **Integration**: Finds TCS-regulated secreted proteins

## Requirements

- Python 3.7+
- pandas
- numpy
- matplotlib
- seaborn
- requests

See `requirements.txt` for full dependencies.

## Future Development

- Integration with real ProteomeXchange datasets (PXD009225, PXD035145, etc.)
- Advanced machine learning models for secretion prediction
- Protein interaction network analysis
- Time-series analysis for infection studies

