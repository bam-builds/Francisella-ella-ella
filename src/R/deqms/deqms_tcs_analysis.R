#!/usr/bin/env Rscript
# ============================================================
# DEqMS Analysis of TCS Proteins in Francisella
# ============================================================
# Differential Expression analysis for proteomics with PSM-dependent
# variance estimation across multiple stress conditions
#
# Author: DEqMS TCS Analysis Pipeline
# Date: 2025
# Dataset: PXD035145 (Francisella stress response proteomics)
# ============================================================

# ============================================================
# LOAD LIBRARIES
# ============================================================
suppressPackageStartupMessages({
  library(DEqMS)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(gridExtra)
})

cat("=============================================================\n")
cat("DEqMS Analysis: Francisella TCS Proteins\n")
cat("=============================================================\n\n")

# ============================================================
# CONFIGURATION
# ============================================================

# TCS protein mapping (locus tag to common name)
TCS_GENES <- c(
  "FTN_1617" = "QseC",     # Sensor kinase (virulence)
  "FTN_1465" = "PmrA",     # Response regulator (LPS modification)
  "FTN_1714" = "KdpD",     # Sensor kinase (osmotic/K+ sensing)
  "FTN_1715" = "KdpE",     # Response regulator (K+ homeostasis)
  "FTN_1452" = "BfpR",     # Response regulator (biofilm)
  "FTN_1453" = "BfpK"      # Sensor kinase (biofilm)
)

# Analysis thresholds
LOG2FC_THRESHOLD <- 1.0     # 2-fold change
PVAL_THRESHOLD <- 0.05      # Adjusted p-value
MIN_VALID_RATIO <- 0.2      # Minimum 20% non-missing values
PSM_MIN <- 2                # Minimum PSM count

# Output directories
OUTPUT_DIRS <- list(
  tables = "results/tables/deqms",
  figures = "results/figures/deqms",
  qc = "results/figures/deqms/qc"
)

# Create output directories
for (dir in OUTPUT_DIRS) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================

#' Log2 transformation with pseudocount
#'
#' @param x Numeric vector or matrix
#' @param pseudocount Value to add before log transformation (default: 1)
#' @return Log2 transformed values
log2_transform <- function(x, pseudocount = 1) {
  log2(x + pseudocount)
}

#' Median normalization
#'
#' @param mat Matrix with samples as columns
#' @return Normalized matrix
median_normalize <- function(mat) {
  medians <- apply(mat, 2, median, na.rm = TRUE)
  global_median <- median(medians)
  normalized <- sweep(mat, 2, medians - global_median, "-")
  return(normalized)
}

#' Filter proteins by missing value threshold
#'
#' @param mat Matrix to filter
#' @param min_ratio Minimum ratio of non-missing values required
#' @return Filtered matrix
filter_missing <- function(mat, min_ratio = 0.2) {
  min_valid <- min_ratio * ncol(mat)
  keep <- rowSums(!is.na(mat)) >= min_valid
  mat[keep, , drop = FALSE]
}

#' Extract significant results
#'
#' @param results DEqMS results data frame
#' @param logfc_threshold Log2 fold change threshold
#' @param pval_threshold Adjusted p-value threshold
#' @return Filtered significant results
extract_significant <- function(results, logfc_threshold = 1.0, pval_threshold = 0.05) {
  results %>%
    filter(
      abs(logFC) > logfc_threshold,
      sca.adj.pval < pval_threshold
    )
}

#' Create volcano plot
#'
#' @param results Results data frame
#' @param title Plot title
#' @param highlight_genes Genes to highlight (optional)
#' @return ggplot object
plot_volcano <- function(results, title, highlight_genes = NULL) {

  # Add significance category
  results$category <- "Not significant"
  results$category[results$logFC > LOG2FC_THRESHOLD & results$sca.adj.pval < PVAL_THRESHOLD] <- "Upregulated"
  results$category[results$logFC < -LOG2FC_THRESHOLD & results$sca.adj.pval < PVAL_THRESHOLD] <- "Downregulated"
  results$category <- factor(results$category,
                             levels = c("Not significant", "Upregulated", "Downregulated"))

  # Base plot
  p <- ggplot(results, aes(x = logFC, y = -log10(sca.P.Value))) +
    geom_point(aes(color = category), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Not significant" = "grey60",
                                  "Upregulated" = "firebrick",
                                  "Downregulated" = "steelblue")) +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
               linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(PVAL_THRESHOLD),
               linetype = "dashed", color = "grey30") +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 P-value",
      color = "Category"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Add gene labels if provided
  if (!is.null(highlight_genes) && "gene" %in% colnames(results)) {
    highlight_data <- results %>%
      filter(gene %in% highlight_genes)

    if (nrow(highlight_data) > 0) {
      p <- p +
        geom_point(data = highlight_data,
                  aes(x = logFC, y = -log10(sca.P.Value)),
                  color = "gold", size = 4, shape = 21, stroke = 1.5) +
        ggrepel::geom_text_repel(
          data = highlight_data,
          aes(label = gene),
          size = 3.5,
          fontface = "bold",
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "grey50"
        )
    }
  }

  return(p)
}

# ============================================================
# STEP 1: DATA IMPORT AND PREPROCESSING
# ============================================================

cat("Step 1: Data Import and Preprocessing\n")
cat("-------------------------------------------------------------\n")

# Check for input files
protein_file <- "data/processed/PXD035145_proteinGroups.csv"
psm_file <- "data/processed/PXD035145_PSM_counts.csv"

if (!file.exists(protein_file)) {
  cat("WARNING: Protein quantification file not found.\n")
  cat("Expected: ", protein_file, "\n")
  cat("Creating example data structure...\n\n")

  # Create example data for demonstration
  set.seed(42)
  n_proteins <- 100
  n_samples <- 24  # 8 conditions × 3 replicates

  # Generate example protein matrix
  protein_raw <- matrix(
    rlnorm(n_proteins * n_samples, meanlog = 20, sdlog = 2),
    nrow = n_proteins,
    dimnames = list(
      c(paste0("FTN_", sprintf("%04d", 1:94)), names(TCS_GENES)),  # Include TCS genes
      paste0("Sample_", 1:n_samples)
    )
  )

  # Add some missing values
  protein_raw[sample(length(protein_raw), length(protein_raw) * 0.1)] <- NA

  # Generate example PSM counts
  psm_counts <- matrix(
    sample(2:50, n_proteins, replace = TRUE),
    nrow = n_proteins,
    dimnames = list(rownames(protein_raw), "PSM_count")
  )

  # Define conditions
  conditions <- factor(rep(c("Control", "Oxidative", "pH_stress", "Heat_42C",
                            "Cold_25C", "Iron_limit", "Nutrient_starv", "Osmotic"),
                          each = 3))

  cat("  Generated example dataset:\n")
  cat("    Proteins: ", n_proteins, "\n")
  cat("    Samples: ", n_samples, "\n")
  cat("    Conditions: ", length(unique(conditions)), "\n\n")

} else {
  # Load real data
  protein_raw <- read.csv(protein_file, row.names = 1)
  psm_counts <- read.csv(psm_file, row.names = 1)

  # Define conditions based on column names
  # This should be adjusted based on actual experimental design
  conditions <- factor(gsub("_[0-9]+$", "", colnames(protein_raw)))

  cat("  Loaded protein data: ", nrow(protein_raw), " proteins × ",
      ncol(protein_raw), " samples\n")
  cat("  Loaded PSM counts for ", nrow(psm_counts), " proteins\n\n")
}

# ============================================================
# STEP 2: DATA TRANSFORMATION AND FILTERING
# ============================================================

cat("Step 2: Data Transformation and Filtering\n")
cat("-------------------------------------------------------------\n")

# Log2 transformation
protein_log2 <- log2_transform(protein_raw)
cat("  Log2 transformation complete\n")

# Filter by missing values
protein_filtered <- filter_missing(protein_log2, min_ratio = MIN_VALID_RATIO)
cat("  Filtered proteins (≥", MIN_VALID_RATIO * 100, "% valid values): ",
    nrow(protein_filtered), "\n")

# Median normalization
protein_norm <- median_normalize(protein_filtered)
cat("  Median normalization complete\n")

# Match PSM counts to filtered proteins
psm_matched <- psm_counts[rownames(protein_norm), , drop = FALSE]

# Handle any remaining mismatches
if (any(is.na(psm_matched))) {
  cat("  WARNING: Some proteins missing PSM counts, using minimum value\n")
  psm_matched[is.na(psm_matched)] <- PSM_MIN
}

cat("  Data preprocessing complete\n\n")

# ============================================================
# STEP 3: QUALITY CONTROL PLOTS
# ============================================================

cat("Step 3: Quality Control\n")
cat("-------------------------------------------------------------\n")

# Sample distribution boxplot
pdf(file.path(OUTPUT_DIRS$qc, "sample_distributions.pdf"), width = 12, height = 6)
par(mar = c(8, 4, 2, 1))
boxplot(protein_norm, las = 2, col = as.numeric(conditions),
        main = "Normalized Protein Intensities",
        ylab = "Log2 Intensity", xlab = "")
legend("topright", legend = levels(conditions),
       fill = 1:length(levels(conditions)), cex = 0.7)
dev.off()

# Missing value heatmap
pdf(file.path(OUTPUT_DIRS$qc, "missing_values.pdf"), width = 10, height = 8)
missing_matrix <- is.na(protein_norm)
pheatmap(
  missing_matrix[1:min(50, nrow(missing_matrix)), ],
  color = c("lightblue", "darkred"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Missing Value Pattern (Top 50 Proteins)",
  legend_labels = c("Observed", "Missing"),
  show_rownames = FALSE
)
dev.off()

cat("  QC plots saved to:", OUTPUT_DIRS$qc, "\n\n")

# ============================================================
# STEP 4: DESIGN MATRIX AND CONTRASTS
# ============================================================

cat("Step 4: Design Matrix and Contrasts\n")
cat("-------------------------------------------------------------\n")

# Create design matrix
design <- model.matrix(~0 + conditions)
colnames(design) <- levels(conditions)

cat("  Design matrix created:\n")
print(table(conditions))
cat("\n")

# Define contrasts (all conditions vs Control)
if ("Control" %in% levels(conditions)) {
  contrast_names <- setdiff(levels(conditions), "Control")
  contrast_formulas <- paste0(contrast_names, " - Control")

  contrasts_matrix <- makeContrasts(
    contrasts = contrast_formulas,
    levels = design
  )

  cat("  Contrasts defined:", ncol(contrasts_matrix), "\n")
  cat("    ", paste(colnames(contrasts_matrix), collapse = "\n     "), "\n\n")
} else {
  stop("Control condition not found in experimental design")
}

# ============================================================
# STEP 5: LINEAR MODELING WITH LIMMA
# ============================================================

cat("Step 5: Linear Modeling (limma)\n")
cat("-------------------------------------------------------------\n")

# Fit linear model
fit1 <- lmFit(protein_norm, design)
cat("  Linear model fitted\n")

# Fit contrasts
fit2 <- contrasts.fit(fit1, contrasts_matrix)
cat("  Contrasts applied\n")

# Empirical Bayes moderation
fit3 <- eBayes(fit2)
cat("  Empirical Bayes moderation complete\n\n")

# ============================================================
# STEP 6: DEqMS VARIANCE CORRECTION
# ============================================================

cat("Step 6: DEqMS Variance Correction\n")
cat("-------------------------------------------------------------\n")

# Add PSM counts to fit object
fit3$count <- psm_matched[, 1]

# Apply DEqMS
fit_deqms <- spectraCounteBayes(fit3)
cat("  DEqMS variance correction applied\n")

# Generate variance-PSM relationship plot
pdf(file.path(OUTPUT_DIRS$qc, "variance_psm_relationship.pdf"),
    width = 10, height = 6)
VarianceBoxplot(fit_deqms, n = 30,
                main = "Variance vs PSM Count Dependency",
                xlab = "PSM Count Bin")
dev.off()

cat("  Variance-PSM plot saved\n\n")

# ============================================================
# STEP 7: EXTRACT RESULTS
# ============================================================

cat("Step 7: Extract Results for All Contrasts\n")
cat("-------------------------------------------------------------\n")

# Extract results for each contrast
results_list <- list()
for (i in 1:ncol(contrasts_matrix)) {
  contrast_name <- colnames(contrasts_matrix)[i]

  # Extract results
  res <- outputResult(fit_deqms, coef_col = i)
  res$gene <- rownames(res)
  res$contrast <- contrast_name

  # Store in list
  results_list[[contrast_name]] <- res

  # Count significant proteins
  n_sig <- sum(res$sca.adj.pval < PVAL_THRESHOLD & abs(res$logFC) > LOG2FC_THRESHOLD,
               na.rm = TRUE)
  cat("  ", contrast_name, ": ", n_sig, " significant proteins\n")
}

# Combine all results
all_results <- bind_rows(results_list)
cat("\n  Total results extracted: ", nrow(all_results), "\n\n")

# Save complete results
write.csv(all_results,
          file.path(OUTPUT_DIRS$tables, "deqms_all_results.csv"),
          row.names = FALSE)
cat("  Saved: deqms_all_results.csv\n\n")

# ============================================================
# STEP 8: TCS PROTEIN-SPECIFIC ANALYSIS
# ============================================================

cat("Step 8: TCS Protein-Specific Analysis\n")
cat("-------------------------------------------------------------\n")

# Extract TCS results
tcs_loci <- names(TCS_GENES)
tcs_results <- all_results %>%
  filter(gene %in% tcs_loci) %>%
  mutate(
    TCS_name = TCS_GENES[gene],
    TCS_type = case_when(
      TCS_name %in% c("QseC", "KdpD", "BfpK") ~ "Sensor Kinase",
      TCS_name %in% c("PmrA", "KdpE", "BfpR") ~ "Response Regulator",
      TRUE ~ "Unknown"
    )
  )

cat("  TCS proteins found: ", length(unique(tcs_results$gene)), "\n")
cat("    ", paste(unique(tcs_results$TCS_name), collapse = ", "), "\n\n")

# Filter for significance
tcs_significant <- tcs_results %>%
  extract_significant(LOG2FC_THRESHOLD, PVAL_THRESHOLD)

cat("  Significant TCS changes: ", nrow(tcs_significant), "\n")

if (nrow(tcs_significant) > 0) {
  cat("\n  Significant TCS protein changes:\n")
  for (i in 1:nrow(tcs_significant)) {
    row <- tcs_significant[i, ]
    direction <- ifelse(row$logFC > 0, "UP", "DOWN")
    cat(sprintf("    %s (%s) in %s: %s %.2f-fold (p=%.2e)\n",
                row$TCS_name, row$TCS_type, row$contrast,
                direction, 2^abs(row$logFC), row$sca.adj.pval))
  }
}

# Save TCS results
write.csv(tcs_results,
          file.path(OUTPUT_DIRS$tables, "TCS_deqms_all_results.csv"),
          row.names = FALSE)

write.csv(tcs_significant,
          file.path(OUTPUT_DIRS$tables, "TCS_deqms_significant.csv"),
          row.names = FALSE)

cat("\n  Saved TCS results to:", OUTPUT_DIRS$tables, "\n\n")

# ============================================================
# STEP 9: VISUALIZATION - VOLCANO PLOTS
# ============================================================

cat("Step 9: Generating Volcano Plots\n")
cat("-------------------------------------------------------------\n")

for (contrast in unique(tcs_results$contrast)) {
  # Get all proteins for this contrast
  contrast_data <- all_results %>%
    filter(contrast == !!contrast)

  # Create volcano plot highlighting TCS proteins
  tcs_genes_in_contrast <- tcs_results %>%
    filter(contrast == !!contrast) %>%
    pull(gene)

  p <- plot_volcano(
    contrast_data,
    title = paste("TCS Proteins:", gsub("_", " ", contrast)),
    highlight_genes = tcs_genes_in_contrast
  )

  # Save plot
  filename <- file.path(OUTPUT_DIRS$figures,
                       paste0("volcano_", contrast, ".pdf"))
  ggsave(filename, p, width = 10, height = 8)

  cat("  Saved:", basename(filename), "\n")
}

cat("\n")

# ============================================================
# STEP 10: HEATMAP OF TCS EXPRESSION PROFILES
# ============================================================

cat("Step 10: TCS Expression Profile Heatmap\n")
cat("-------------------------------------------------------------\n")

# Extract normalized expression for TCS proteins
tcs_in_data <- intersect(tcs_loci, rownames(protein_norm))

if (length(tcs_in_data) > 0) {
  tcs_expression <- protein_norm[tcs_in_data, , drop = FALSE]
  rownames(tcs_expression) <- TCS_GENES[tcs_in_data]

  # Center and scale by row
  tcs_centered <- t(scale(t(tcs_expression)))

  # Create annotation
  annotation_col <- data.frame(
    Condition = conditions,
    row.names = colnames(tcs_centered)
  )

  # Color palette
  condition_colors <- brewer.pal(
    min(length(unique(conditions)), 9),
    "Set1"
  )
  names(condition_colors) <- levels(conditions)

  annotation_colors <- list(
    Condition = condition_colors
  )

  # Generate heatmap
  pdf(file.path(OUTPUT_DIRS$figures, "TCS_expression_heatmap.pdf"),
      width = 12, height = 8)

  pheatmap(
    tcs_centered,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    clustering_distance_rows = "correlation",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-2, 2, length.out = 101),
    main = "TCS Protein Expression Across Stress Conditions",
    fontsize = 10,
    fontsize_row = 11,
    fontsize_col = 9
  )

  dev.off()

  cat("  Heatmap saved: TCS_expression_heatmap.pdf\n\n")
} else {
  cat("  WARNING: No TCS proteins found in normalized data\n\n")
}

# ============================================================
# STEP 11: SUMMARY STATISTICS AND REPORT
# ============================================================

cat("Step 11: Summary Report\n")
cat("=============================================================\n\n")

# Overall statistics
summary_stats <- list(
  total_proteins = nrow(protein_norm),
  total_samples = ncol(protein_norm),
  conditions_tested = length(unique(conditions)),
  contrasts_analyzed = ncol(contrasts_matrix),
  tcs_proteins_detected = length(tcs_in_data),
  tcs_significant_changes = nrow(tcs_significant)
)

cat("ANALYSIS SUMMARY\n")
cat("-------------------------------------------------------------\n")
cat(sprintf("Total proteins analyzed:            %d\n", summary_stats$total_proteins))
cat(sprintf("Total samples:                      %d\n", summary_stats$total_samples))
cat(sprintf("Experimental conditions:            %d\n", summary_stats$conditions_tested))
cat(sprintf("Contrasts analyzed:                 %d\n", summary_stats$contrasts_analyzed))
cat(sprintf("TCS proteins detected:              %d / %d\n",
            summary_stats$tcs_proteins_detected, length(TCS_GENES)))
cat(sprintf("Significant TCS changes:            %d\n", summary_stats$tcs_significant_changes))
cat("\n")

# Per-contrast summary
cat("DIFFERENTIAL EXPRESSION BY CONTRAST\n")
cat("-------------------------------------------------------------\n")

contrast_summary <- all_results %>%
  group_by(contrast) %>%
  summarise(
    total = n(),
    sig_up = sum(logFC > LOG2FC_THRESHOLD & sca.adj.pval < PVAL_THRESHOLD, na.rm = TRUE),
    sig_down = sum(logFC < -LOG2FC_THRESHOLD & sca.adj.pval < PVAL_THRESHOLD, na.rm = TRUE),
    .groups = "drop"
  )

print(contrast_summary)
cat("\n")

# Save summary
write.csv(contrast_summary,
          file.path(OUTPUT_DIRS$tables, "deqms_summary_statistics.csv"),
          row.names = FALSE)

# ============================================================
# COMPLETION
# ============================================================

cat("=============================================================\n")
cat("DEqMS ANALYSIS COMPLETE\n")
cat("=============================================================\n\n")

cat("Output files generated:\n")
cat("  Tables:  ", OUTPUT_DIRS$tables, "\n")
cat("  Figures: ", OUTPUT_DIRS$figures, "\n")
cat("  QC:      ", OUTPUT_DIRS$qc, "\n\n")

cat("Key results files:\n")
cat("  - deqms_all_results.csv\n")
cat("  - TCS_deqms_significant.csv\n")
cat("  - TCS_expression_heatmap.pdf\n")
cat("  - Volcano plots for each contrast\n\n")

cat("Analysis completed successfully!\n")
cat("=============================================================\n")
