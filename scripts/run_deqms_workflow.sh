#!/bin/bash
# ============================================================
# Master Script for DEqMS TCS Analysis Workflow
# ============================================================
# Runs the complete analysis pipeline from data preparation
# to final visualization
#
# Usage: bash scripts/run_deqms_workflow.sh
# ============================================================

set -e  # Exit on error

echo "============================================================="
echo "DEqMS TCS Analysis - Complete Workflow"
echo "============================================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "Working directory: $PROJECT_DIR"
echo ""

# ============================================================
# Step 0: Check R installation
# ============================================================

echo "Step 0: Checking R installation"
echo "-------------------------------------------------------------"

if ! command -v Rscript &> /dev/null; then
    echo -e "${RED}ERROR: R is not installed or not in PATH${NC}"
    echo "Please install R (version >= 4.0.0) and try again"
    exit 1
fi

R_VERSION=$(Rscript --version 2>&1 | grep -oP 'version \K[0-9.]+' | head -1)
echo "R version: $R_VERSION"
echo ""

# ============================================================
# Step 1: Install R packages (optional)
# ============================================================

if [ "$1" == "--install-packages" ]; then
    echo "Step 1: Installing R packages"
    echo "-------------------------------------------------------------"
    Rscript scripts/install_r_packages.R

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Package installation successful${NC}"
    else
        echo -e "${RED}✗ Package installation failed${NC}"
        exit 1
    fi
    echo ""
fi

# ============================================================
# Step 2: Check for input data
# ============================================================

echo "Step 2: Checking input data"
echo "-------------------------------------------------------------"

DATA_EXISTS=false

# Check for MaxQuant data
if [ -f "data/raw/maxquant/proteinGroups.txt" ]; then
    echo "✓ Found MaxQuant proteinGroups.txt"
    DATA_EXISTS=true
elif [ -f "data/processed/PXD035145_proteinGroups.csv" ]; then
    echo "✓ Found processed protein data"
    DATA_EXISTS=true
else
    echo -e "${YELLOW}⚠ No input data found${NC}"
    echo "The analysis will run with example data"
fi

echo ""

# ============================================================
# Step 3: Prepare data (if MaxQuant files exist)
# ============================================================

if [ -f "data/raw/maxquant/proteinGroups.txt" ]; then
    echo "Step 3: Preparing MaxQuant data"
    echo "-------------------------------------------------------------"

    Rscript src/R/deqms/prepare_maxquant_data.R

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Data preparation successful${NC}"
    else
        echo -e "${RED}✗ Data preparation failed${NC}"
        exit 1
    fi
    echo ""
else
    echo "Step 3: Skipping data preparation (no MaxQuant files)"
    echo "-------------------------------------------------------------"
    echo "Using example data or existing processed files"
    echo ""
fi

# ============================================================
# Step 4: Run DEqMS analysis
# ============================================================

echo "Step 4: Running DEqMS analysis"
echo "-------------------------------------------------------------"

Rscript src/R/deqms/deqms_tcs_analysis.R

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ DEqMS analysis completed successfully${NC}"
else
    echo -e "${RED}✗ DEqMS analysis failed${NC}"
    exit 1
fi
echo ""

# ============================================================
# Step 5: Generate visualizations
# ============================================================

echo "Step 5: Generating publication figures"
echo "-------------------------------------------------------------"

Rscript src/R/deqms/visualize_tcs_results.R

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Visualization completed successfully${NC}"
else
    echo -e "${YELLOW}⚠ Visualization failed (analysis results still available)${NC}"
fi
echo ""

# ============================================================
# Step 6: Summary
# ============================================================

echo "============================================================="
echo "WORKFLOW COMPLETE"
echo "============================================================="
echo ""

echo "Output directories:"
echo "  📊 Tables:  results/tables/deqms/"
echo "  📈 Figures: results/figures/deqms/"
echo "  🔍 QC:      results/figures/deqms/qc/"
echo "  📰 Publication: results/figures/deqms/publication/"
echo ""

echo "Key result files:"
if [ -f "results/tables/deqms/TCS_deqms_significant.csv" ]; then
    N_SIG=$(tail -n +2 "results/tables/deqms/TCS_deqms_significant.csv" | wc -l)
    echo "  ✓ TCS_deqms_significant.csv ($N_SIG significant changes)"
else
    echo "  ⚠ TCS_deqms_significant.csv (not found)"
fi

if [ -f "results/figures/deqms/TCS_expression_heatmap.pdf" ]; then
    echo "  ✓ TCS_expression_heatmap.pdf"
fi

if [ -f "results/figures/deqms/publication/Fig1_TCS_heatmap_logFC.pdf" ]; then
    echo "  ✓ Publication-ready figures (6 figures generated)"
fi

echo ""
echo "Next steps:"
echo "  1. Review QC plots in results/figures/deqms/qc/"
echo "  2. Check significant TCS changes in TCS_deqms_significant.csv"
echo "  3. Examine publication figures for biological insights"
echo "  4. See docs/DEqMS_TCS_Analysis_Guide.md for interpretation"
echo ""

echo "============================================================="
echo "Analysis completed successfully!"
echo "============================================================="
