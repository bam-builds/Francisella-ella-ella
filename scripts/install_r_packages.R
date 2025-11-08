#!/usr/bin/env Rscript
# ============================================================
# R Package Installation for DEqMS TCS Analysis
# ============================================================
# Installs all required R packages for the DEqMS workflow
#
# Usage: Rscript scripts/install_r_packages.R
# ============================================================

cat("=============================================================\n")
cat("Installing R Packages for DEqMS TCS Analysis\n")
cat("=============================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ============================================================
# PACKAGE LISTS
# ============================================================

# CRAN packages
cran_packages <- c(
  "dplyr",           # Data manipulation
  "tidyr",           # Data tidying
  "ggplot2",         # Plotting
  "pheatmap",        # Heatmaps
  "RColorBrewer",    # Color palettes
  "gridExtra",       # Multiple plots
  "ggrepel",         # Repel overlapping text labels
  "stringr",         # String manipulation
  "readr",           # Fast file reading
  "tibble"           # Modern data frames
)

# Bioconductor packages
bioc_packages <- c(
  "limma",           # Linear models for microarray data
  "DEqMS",           # Differential expression for proteomics
  "ComplexHeatmap",  # Advanced heatmaps
  "circlize"         # Circular visualization
)

# ============================================================
# HELPER FUNCTION: CHECK AND INSTALL
# ============================================================

#' Check if package is installed, install if not
#'
#' @param pkg Package name
#' @param source Source ("CRAN" or "Bioconductor")
#' @return TRUE if successful
install_if_missing <- function(pkg, source = "CRAN") {

  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  ✓", pkg, "- already installed\n")
    return(TRUE)
  }

  cat("  Installing", pkg, "from", source, "...\n")

  tryCatch({
    if (source == "CRAN") {
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", quiet = TRUE)
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    }

    # Verify installation
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("  ✓", pkg, "- installed successfully\n")
      return(TRUE)
    } else {
      cat("  ✗", pkg, "- installation failed\n")
      return(FALSE)
    }

  }, error = function(e) {
    cat("  ✗", pkg, "- ERROR:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# ============================================================
# INSTALL CRAN PACKAGES
# ============================================================

cat("Step 1: Installing CRAN Packages\n")
cat("-------------------------------------------------------------\n")

cran_success <- sapply(cran_packages, function(pkg) {
  install_if_missing(pkg, source = "CRAN")
})

cat("\nCRAN packages installed:", sum(cran_success), "/", length(cran_packages), "\n\n")

# ============================================================
# INSTALL BIOCONDUCTOR PACKAGES
# ============================================================

cat("Step 2: Installing Bioconductor Packages\n")
cat("-------------------------------------------------------------\n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("  Installing BiocManager...\n")
  install.packages("BiocManager", quiet = TRUE)
}

bioc_success <- sapply(bioc_packages, function(pkg) {
  install_if_missing(pkg, source = "Bioconductor")
})

cat("\nBioconductor packages installed:", sum(bioc_success), "/", length(bioc_packages), "\n\n")

# ============================================================
# VERIFICATION
# ============================================================

cat("Step 3: Verification\n")
cat("=============================================================\n\n")

all_packages <- c(cran_packages, bioc_packages)
installed <- sapply(all_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

cat("Package Installation Summary:\n")
cat("-------------------------------------------------------------\n")
cat("Total packages required: ", length(all_packages), "\n")
cat("Successfully installed:  ", sum(installed), "\n")
cat("Failed installations:    ", sum(!installed), "\n\n")

if (sum(!installed) > 0) {
  cat("Failed packages:\n")
  cat("  ", paste(names(installed)[!installed], collapse = "\n   "), "\n\n")
  cat("NOTE: Please install these packages manually\n\n")
}

# ============================================================
# TEST CRITICAL PACKAGES
# ============================================================

cat("Step 4: Testing Critical Packages\n")
cat("-------------------------------------------------------------\n")

critical_packages <- c("DEqMS", "limma", "ggplot2", "dplyr")

all_critical_ok <- TRUE

for (pkg in critical_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    # Try to load the package
    tryCatch({
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat("  ✓", pkg, "- loaded successfully\n")
    }, error = function(e) {
      cat("  ✗", pkg, "- failed to load:", conditionMessage(e), "\n")
      all_critical_ok <- FALSE
    })
  } else {
    cat("  ✗", pkg, "- not installed\n")
    all_critical_ok <- FALSE
  }
}

cat("\n")

# ============================================================
# FINAL STATUS
# ============================================================

cat("=============================================================\n")

if (all_critical_ok && sum(!installed) == 0) {
  cat("SUCCESS: All packages installed and working!\n")
  cat("=============================================================\n\n")
  cat("You can now run:\n")
  cat("  Rscript src/R/deqms/deqms_tcs_analysis.R\n\n")
  quit(save = "no", status = 0)

} else if (all_critical_ok) {
  cat("PARTIAL SUCCESS: Critical packages installed\n")
  cat("=============================================================\n\n")
  cat("Some optional packages failed, but you can proceed.\n")
  cat("Missing packages may affect visualization features.\n\n")
  quit(save = "no", status = 0)

} else {
  cat("INSTALLATION INCOMPLETE\n")
  cat("=============================================================\n\n")
  cat("Critical packages are missing. Please install manually:\n\n")
  for (pkg in critical_packages[!critical_packages %in% names(installed)[installed]]) {
    cat("  install.packages('", pkg, "')\n", sep = "")
  }
  cat("\nFor Bioconductor packages:\n")
  cat("  BiocManager::install(c('DEqMS', 'limma'))\n\n")
  quit(save = "no", status = 1)
}
