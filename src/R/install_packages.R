# ============================================================================
# R Package Installation Script for TCS Proteomics Visualization
# ============================================================================

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Installing %s...\n", pkg))
      install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    } else {
      cat(sprintf("✓ %s already installed\n", pkg))
    }
  }
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(packages) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }

  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Installing %s from Bioconductor...\n", pkg))
      BiocManager::install(pkg, update = FALSE)
    } else {
      cat(sprintf("✓ %s already installed\n", pkg))
    }
  }
}

# ============================================================================
# Core Data Manipulation and Visualization Packages
# ============================================================================

cran_packages <- c(
  # Core tidyverse
  "tidyverse",      # Data manipulation and visualization
  "dplyr",          # Data manipulation
  "ggplot2",        # Plotting
  "tidyr",          # Data tidying

  # Visualization packages
  "viridis",        # Color scales
  "RColorBrewer",   # Color palettes
  "ggrepel",        # Label repulsion in plots
  "cowplot",        # Plot composition
  "ggsci",          # Scientific journal color palettes
  "scales",         # Scale functions for visualization
  "gridExtra",      # Grid arrangements

  # Network visualization
  "network",        # Network objects
  "ggnetwork",      # Network visualization with ggplot2

  # Heatmap packages
  "pheatmap",       # Pretty heatmaps
  "circlize"        # Circular visualization and color mapping
)

# ============================================================================
# Bioconductor Packages
# ============================================================================

bioc_packages <- c(
  "ComplexHeatmap"  # Advanced heatmap visualization
)

# ============================================================================
# Installation
# ============================================================================

cat("================================================================================\n")
cat("Installing R Packages for TCS Proteomics Visualization\n")
cat("================================================================================\n\n")

cat("Step 1: Installing CRAN packages...\n")
cat("--------------------------------------------------------------------------------\n")
install_if_missing(cran_packages)

cat("\nStep 2: Installing Bioconductor packages...\n")
cat("--------------------------------------------------------------------------------\n")
install_bioc_if_missing(bioc_packages)

# ============================================================================
# Verification
# ============================================================================

cat("\n================================================================================\n")
cat("Package Installation Complete!\n")
cat("================================================================================\n\n")

cat("Verifying installations...\n")
all_packages <- c(cran_packages, bioc_packages)
missing_packages <- c()

for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  cat("\n✓ All packages successfully installed!\n")
  cat("\nYou can now run the visualization script:\n")
  cat("  Rscript src/R/tcs_proteomics_visualization.R\n\n")
} else {
  cat("\n⚠ The following packages failed to install:\n")
  for (pkg in missing_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease try installing these packages manually.\n\n")
}

# ============================================================================
# Session Info
# ============================================================================

cat("R Session Information:\n")
cat("--------------------------------------------------------------------------------\n")
print(sessionInfo())
