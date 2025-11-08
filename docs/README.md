# Francisella TCS Analysis Documentation

This directory contains comprehensive documentation for the Francisella Two-Component System (TCS) proteomics analysis pipeline.

## Available Documentation

### [DEqMS TCS Analysis Guide](DEqMS_TCS_Analysis_Guide.md)

**Complete guide for DEqMS-based differential expression analysis**

Topics covered:
- Installation and setup
- Data preparation from MaxQuant output
- Running the analysis workflow
- Interpreting results
- Biological significance of TCS proteins
- Troubleshooting common issues
- Advanced customization

**When to use:** For quantitative proteomics analysis of TCS protein expression across multiple stress conditions using DEqMS methodology.

---

## Quick Start

### For DEqMS Analysis

```bash
# 1. Install R packages (one-time setup)
bash scripts/run_deqms_workflow.sh --install-packages

# 2. Place MaxQuant files in data/raw/maxquant/
#    - proteinGroups.txt
#    - evidence.txt

# 3. Run complete workflow
bash scripts/run_deqms_workflow.sh
```

### For Python-based Analysis

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Run analysis
python analysis.py
```

---

## Analysis Workflows

### 1. DEqMS Differential Expression (R)

**Purpose:** Statistical analysis of TCS protein abundance changes

**Input:**
- MaxQuant proteinGroups.txt
- MaxQuant evidence.txt

**Output:**
- Differential expression statistics
- Volcano plots
- Expression heatmaps
- QC reports

**Scripts:**
- `src/R/deqms/deqms_tcs_analysis.R`
- `src/R/deqms/prepare_maxquant_data.R`
- `src/R/deqms/visualize_tcs_results.R`

### 2. Python Proteomics Pipeline

**Purpose:** Comprehensive TCS analysis including PTMs and secretion

**Input:**
- Proteomics data files
- Phosphoproteomics data
- Protein sequences

**Output:**
- TCS protein identification
- Phosphorylation site analysis
- Secretion signal prediction
- Integrated analysis

**Scripts:**
- `analysis.py`
- `src/python/data_retrieval.py`
- `src/python/ptm_analysis.py`
- `src/python/secretion_prediction.py`

---

## TCS Proteins Analyzed

### QseC/PmrA System
- **Function:** Master virulence regulator, controls FPI
- **Genes:** QseC (FTN_1617), PmrA (FTN_1465)
- **Expected activation:** Oxidative stress, iron limitation, nutrient starvation

### KdpD/KdpE System
- **Function:** Potassium homeostasis
- **Genes:** KdpD (FTN_1714), KdpE (FTN_1715)
- **Expected activation:** Osmotic stress, low K+
- **Note:** Pseudogene in SchuS4, functional in F. novicida

### BfpR/BfpK System
- **Function:** Biofilm regulation
- **Genes:** BfpK (FTN_1453), BfpR (FTN_1452)
- **Expected activation:** Mg2+ limitation, antimicrobial peptides

---

## Data Sources

### Primary Dataset
- **Accession:** PXD035145
- **Title:** Francisella stress response proteomics
- **Description:** Multi-condition proteomics analysis of Francisella novicida

### Additional Resources
- ProteomeXchange repository
- UniProt (Francisella proteome)
- KEGG pathways (TCS systems)

---

## File Organization

```
docs/
├── README.md                      # This file
└── DEqMS_TCS_Analysis_Guide.md   # Complete DEqMS guide
```

---

## Getting Help

### For DEqMS-specific questions
- See [DEqMS_TCS_Analysis_Guide.md](DEqMS_TCS_Analysis_Guide.md)
- Troubleshooting section covers common issues
- References to original DEqMS publication

### For general proteomics questions
- Check MaxQuant documentation
- Review limma user guide for statistical methods

### For Francisella biology
- References in main documentation
- TCS protein function descriptions in analysis guide

---

## Citation

If you use this pipeline, please cite the appropriate methods:

**For DEqMS analysis:**
> Zhu Y, et al. DEqMS: A Method for Accurate Variance Estimation in Differential Protein Expression Analysis. Mol Cell Proteomics. 2020;19(6):1047-1057.

**For limma:**
> Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47.

---

## Updates

- **2025:** Initial DEqMS workflow implementation
- **2025:** Added comprehensive TCS analysis documentation
- **2025:** Multi-stress condition comparison pipeline

---

## License

This documentation and associated code are provided for academic research purposes.
