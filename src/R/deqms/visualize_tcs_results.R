#!/usr/bin/env Rscript
# ============================================================
# Advanced Visualization for TCS DEqMS Results
# ============================================================
# Creates publication-quality figures for TCS protein
# differential expression analysis
#
# Author: DEqMS TCS Analysis Pipeline
# Date: 2025
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(ggrepel)
  library(ComplexHeatmap)  # For advanced heatmaps
  library(circlize)        # For color mapping
})

cat("=============================================================\n")
cat("TCS DEqMS Results Visualization\n")
cat("=============================================================\n\n")

# ============================================================
# CONFIGURATION
# ============================================================

RESULTS_DIR <- "results/tables/deqms"
OUTPUT_DIR <- "results/figures/deqms/publication"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load results
TCS_RESULTS_FILE <- file.path(RESULTS_DIR, "TCS_deqms_all_results.csv")
TCS_SIG_FILE <- file.path(RESULTS_DIR, "TCS_deqms_significant.csv")

if (!file.exists(TCS_RESULTS_FILE)) {
  stop("ERROR: TCS results file not found. Run deqms_tcs_analysis.R first.")
}

tcs_results <- read.csv(TCS_RESULTS_FILE, stringsAsFactors = FALSE)
tcs_sig <- read.csv(TCS_SIG_FILE, stringsAsFactors = FALSE)

cat("Loaded results:\n")
cat("  Total TCS observations: ", nrow(tcs_results), "\n")
cat("  Significant changes:    ", nrow(tcs_sig), "\n\n")

# ============================================================
# FIGURE 1: TCS Response Heatmap Across All Conditions
# ============================================================

cat("Creating Figure 1: TCS Response Heatmap\n")
cat("-------------------------------------------------------------\n")

# Reshape data for heatmap
tcs_matrix <- tcs_results %>%
  select(TCS_name, contrast, logFC) %>%
  spread(key = contrast, value = logFC) %>%
  as.data.frame()

rownames(tcs_matrix) <- tcs_matrix$TCS_name
tcs_matrix <- tcs_matrix[, -1]
tcs_matrix <- as.matrix(tcs_matrix)

# Remove columns/rows with all NA
tcs_matrix <- tcs_matrix[rowSums(!is.na(tcs_matrix)) > 0,
                         colSums(!is.na(tcs_matrix)) > 0,
                         drop = FALSE]

# Create annotation for TCS type
tcs_annotation <- data.frame(
  Type = c(
    "QseC" = "Sensor Kinase",
    "PmrA" = "Response Regulator",
    "KdpD" = "Sensor Kinase",
    "KdpE" = "Response Regulator",
    "BfpR" = "Response Regulator",
    "BfpK" = "Sensor Kinase"
  )[rownames(tcs_matrix)],
  row.names = rownames(tcs_matrix)
)

# Colors
type_colors <- c(
  "Sensor Kinase" = "#E69F00",
  "Response Regulator" = "#56B4E9"
)

annotation_colors <- list(
  Type = type_colors
)

# Create heatmap
pdf(file.path(OUTPUT_DIR, "Fig1_TCS_heatmap_logFC.pdf"),
    width = 10, height = 6)

pheatmap(
  tcs_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  annotation_row = tcs_annotation,
  annotation_colors = annotation_colors,
  main = "TCS Protein Fold Changes Across Stress Conditions",
  fontsize = 11,
  fontsize_row = 12,
  fontsize_col = 10,
  na_col = "grey90",
  border_color = NA,
  cellwidth = 40,
  cellheight = 30
)

dev.off()

cat("  Saved: Fig1_TCS_heatmap_logFC.pdf\n\n")

# ============================================================
# FIGURE 2: Multi-Panel Volcano Plots
# ============================================================

cat("Creating Figure 2: Multi-Panel Volcano Plots\n")
cat("-------------------------------------------------------------\n")

# Get unique contrasts
contrasts <- unique(tcs_results$contrast)

# Create volcano plots for each contrast
volcano_plots <- lapply(contrasts, function(contrast_name) {

  # Get data for this contrast
  contrast_data <- tcs_results %>%
    filter(contrast == contrast_name)

  # Determine significance
  contrast_data$sig <- ifelse(
    abs(contrast_data$logFC) > 1 & contrast_data$sca.adj.pval < 0.05,
    "Significant",
    "Not Significant"
  )

  # Create plot
  p <- ggplot(contrast_data, aes(x = logFC, y = -log10(sca.P.Value))) +
    geom_point(aes(color = sig, size = sig), alpha = 0.7) +
    scale_color_manual(values = c("Not Significant" = "grey70",
                                  "Significant" = "red")) +
    scale_size_manual(values = c("Not Significant" = 2,
                                 "Significant" = 4)) +
    geom_text_repel(
      data = contrast_data %>% filter(sig == "Significant"),
      aes(label = TCS_name),
      size = 3.5,
      fontface = "bold",
      box.padding = 0.5,
      max.overlaps = 20
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(
      title = gsub("_", " ", contrast_name),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    xlim(c(-max(abs(contrast_data$logFC), na.rm = TRUE) - 0.5,
           max(abs(contrast_data$logFC), na.rm = TRUE) + 0.5))

  return(p)
})

# Arrange in grid
pdf(file.path(OUTPUT_DIR, "Fig2_volcano_panels.pdf"),
    width = 16, height = 12)

do.call(grid.arrange, c(volcano_plots, ncol = 3))

dev.off()

cat("  Saved: Fig2_volcano_panels.pdf\n\n")

# ============================================================
# FIGURE 3: TCS Protein-Specific Expression Profiles
# ============================================================

cat("Creating Figure 3: Individual TCS Expression Profiles\n")
cat("-------------------------------------------------------------\n")

# Create line plots for each TCS protein
tcs_profile_plots <- lapply(unique(tcs_results$TCS_name), function(tcs_name) {

  tcs_data <- tcs_results %>%
    filter(TCS_name == tcs_name) %>%
    mutate(
      significant = abs(logFC) > 1 & sca.adj.pval < 0.05,
      contrast_clean = gsub("_", " ", contrast)
    )

  # Determine protein type
  protein_type <- unique(tcs_data$TCS_type)[1]
  color_choice <- ifelse(protein_type == "Sensor Kinase", "#E69F00", "#56B4E9")

  p <- ggplot(tcs_data, aes(x = reorder(contrast_clean, logFC), y = logFC)) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = c(-1, 1), color = "grey60", linetype = "dashed") +
    geom_col(aes(fill = significant), width = 0.7) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = color_choice)) +
    geom_errorbar(
      aes(ymin = logFC - sca.SE, ymax = logFC + sca.SE),
      width = 0.3,
      color = "black"
    ) +
    labs(
      title = paste(tcs_name, "-", protein_type),
      x = NULL,
      y = "Log2 Fold Change"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    ) +
    coord_flip()

  return(p)
})

# Save individual TCS profiles
pdf(file.path(OUTPUT_DIR, "Fig3_TCS_expression_profiles.pdf"),
    width = 10, height = 12)

do.call(grid.arrange, c(tcs_profile_plots, ncol = 2))

dev.off()

cat("  Saved: Fig3_TCS_expression_profiles.pdf\n\n")

# ============================================================
# FIGURE 4: Significance Summary Barplot
# ============================================================

cat("Creating Figure 4: Significance Summary\n")
cat("-------------------------------------------------------------\n")

# Count significant changes per TCS protein
sig_summary <- tcs_results %>%
  mutate(
    direction = case_when(
      logFC > 1 & sca.adj.pval < 0.05 ~ "Upregulated",
      logFC < -1 & sca.adj.pval < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  filter(direction != "Not Significant") %>%
  group_by(TCS_name, TCS_type, direction) %>%
  summarise(count = n(), .groups = "drop")

# Create stacked bar plot
p4 <- ggplot(sig_summary, aes(x = reorder(TCS_name, count), y = count, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Upregulated" = "firebrick",
                               "Downregulated" = "steelblue")) +
  labs(
    title = "Number of Significant Changes per TCS Protein",
    x = "TCS Protein",
    y = "Number of Significant Contrasts",
    fill = "Direction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  coord_flip()

ggsave(file.path(OUTPUT_DIR, "Fig4_significance_summary.pdf"),
       p4, width = 8, height = 6)

cat("  Saved: Fig4_significance_summary.pdf\n\n")

# ============================================================
# FIGURE 5: P-value Distribution
# ============================================================

cat("Creating Figure 5: P-value Distribution\n")
cat("-------------------------------------------------------------\n")

p5 <- ggplot(tcs_results, aes(x = sca.P.Value)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~ TCS_name, scales = "free_y", ncol = 3) +
  labs(
    title = "P-value Distribution by TCS Protein",
    x = "P-value",
    y = "Frequency"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(file.path(OUTPUT_DIR, "Fig5_pvalue_distribution.pdf"),
       p5, width = 12, height = 8)

cat("  Saved: Fig5_pvalue_distribution.pdf\n\n")

# ============================================================
# FIGURE 6: Correlation Between TCS Proteins
# ============================================================

cat("Creating Figure 6: TCS Protein Correlation\n")
cat("-------------------------------------------------------------\n")

# Create correlation matrix of log fold changes
tcs_cor_data <- tcs_results %>%
  select(TCS_name, contrast, logFC) %>%
  spread(key = TCS_name, value = logFC)

rownames(tcs_cor_data) <- tcs_cor_data$contrast
tcs_cor_data <- tcs_cor_data[, -1]

# Calculate correlation
tcs_cor <- cor(tcs_cor_data, use = "pairwise.complete.obs")

# Create correlation heatmap
pdf(file.path(OUTPUT_DIR, "Fig6_TCS_correlation.pdf"),
    width = 8, height = 7)

pheatmap(
  tcs_cor,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Correlation of TCS Protein Responses",
  fontsize = 11,
  fontsize_number = 9,
  border_color = "grey60"
)

dev.off()

cat("  Saved: Fig6_TCS_correlation.pdf\n\n")

# ============================================================
# SUMMARY TABLE FOR PUBLICATION
# ============================================================

cat("Creating Supplementary Table\n")
cat("-------------------------------------------------------------\n")

# Create formatted results table
pub_table <- tcs_results %>%
  select(
    `TCS Protein` = TCS_name,
    `Protein Type` = TCS_type,
    Contrast = contrast,
    `Log2 FC` = logFC,
    `P-value` = sca.P.Value,
    `Adj. P-value` = sca.adj.pval,
    `Fold Change` = logFC
  ) %>%
  mutate(
    `Log2 FC` = round(`Log2 FC`, 2),
    `Fold Change` = round(2^abs(`Fold Change`), 2),
    `P-value` = formatC(`P-value`, format = "e", digits = 2),
    `Adj. P-value` = formatC(`Adj. P-value`, format = "e", digits = 2),
    Contrast = gsub("_", " ", Contrast)
  ) %>%
  arrange(`TCS Protein`, Contrast)

write.csv(pub_table,
          file.path(OUTPUT_DIR, "TableS1_TCS_DEqMS_results.csv"),
          row.names = FALSE)

cat("  Saved: TableS1_TCS_DEqMS_results.csv\n\n")

# ============================================================
# BIOLOGICAL INTERPRETATION SUMMARY
# ============================================================

cat("Generating Biological Interpretation Summary\n")
cat("-------------------------------------------------------------\n")

# Expected biological patterns
bio_summary <- list(
  QseC_PmrA = c("Oxidative", "Iron_limit", "Nutrient_starv"),
  KdpD_KdpE = c("Osmotic", "Cold_25C"),
  BfpR_BfpK = c("Osmotic")
)

cat("\nBiological Pattern Analysis:\n\n")

for (system in names(bio_summary)) {
  cat("System:", system, "\n")

  expected_conditions <- bio_summary[[system]]
  proteins <- strsplit(system, "_")[[1]]

  for (protein in proteins) {
    protein_data <- tcs_sig %>%
      filter(TCS_name == protein) %>%
      mutate(contrast_short = gsub("_.*", "", contrast))

    if (nrow(protein_data) > 0) {
      cat("  ", protein, ":\n")
      for (condition in expected_conditions) {
        matches <- protein_data %>%
          filter(grepl(condition, contrast, ignore.case = TRUE))

        if (nrow(matches) > 0) {
          cat("    ✓", condition, "- SIGNIFICANT (",
              sprintf("%.2f-fold, p=%.2e",
                      2^abs(matches$logFC[1]),
                      matches$sca.adj.pval[1]), ")\n")
        } else {
          cat("    ✗", condition, "- not significant\n")
        }
      }
    }
  }
  cat("\n")
}

# ============================================================
# COMPLETION
# ============================================================

cat("=============================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("=============================================================\n\n")

cat("Generated figures:\n")
cat("  1. TCS_heatmap_logFC.pdf         - Overview heatmap\n")
cat("  2. volcano_panels.pdf            - Multi-panel volcano plots\n")
cat("  3. TCS_expression_profiles.pdf   - Individual TCS profiles\n")
cat("  4. significance_summary.pdf      - Summary barplot\n")
cat("  5. pvalue_distribution.pdf       - P-value distributions\n")
cat("  6. TCS_correlation.pdf           - Correlation heatmap\n\n")

cat("Supplementary materials:\n")
cat("  - TableS1_TCS_DEqMS_results.csv\n\n")

cat("All files saved to:", OUTPUT_DIR, "\n\n")

cat("=============================================================\n")
