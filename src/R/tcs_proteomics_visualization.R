# ============================================================================
# FRANCISELLA TCS PROTEOMICS DATA VISUALIZATION
# Focus: Up/Down Regulation of Two-Component Systems
# ============================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(cowplot)
library(ggsci)
library(scales)
library(network)
library(ggnetwork)

# Set publication-quality theme
theme_set(theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  ))

# ============================================================================
# SECTION 1: MOCK DATA GENERATION (Replace with actual ProteomeXchange data)
# ============================================================================

generate_mock_tcs_data <- function() {

  # Define TCS proteins and related genes
  tcs_proteins <- c("QseB/PmrA", "QseC", "KdpE", "KdpD", "BfpR", "BfpK", "MglA", "SspA")
  fpi_genes <- c("IglA", "IglB", "IglC", "IglD", "VgrG", "IglJ", "PdpA", "PdpB")
  stress_genes <- c("KatG", "AhpC", "SodB", "OxyR", "TrxA", "TrxB", "GrxA", "GshA")

  all_proteins <- c(tcs_proteins, fpi_genes, stress_genes)

  # Generate expression data for different conditions
  conditions <- list(
    "Control" = rnorm(length(all_proteins), 0, 0.5),
    "Oxidative_2h" = c(rnorm(8, 1.5, 0.5), rnorm(8, 2.5, 0.5), rnorm(8, 3, 0.5)),
    "Oxidative_8h" = c(rnorm(8, -1, 0.5), rnorm(8, 3.5, 0.5), rnorm(8, 2, 0.5)),
    "pH5.5" = c(rnorm(8, 2, 0.5), rnorm(8, -1, 0.5), rnorm(8, 1, 0.5)),
    "pH8.5" = c(rnorm(8, -1.5, 0.5), rnorm(8, 0.5, 0.5), rnorm(8, -0.5, 0.5)),
    "Temp_37C" = c(rnorm(8, 1, 0.5), rnorm(8, 2, 0.5), rnorm(8, 0.5, 0.5)),
    "Iron_limit" = c(rnorm(8, -2, 0.5), rnorm(8, 1, 0.5), rnorm(8, 2.5, 0.5)),
    "Macrophage_2h" = c(rnorm(8, 3, 0.5), rnorm(8, 4, 0.5), rnorm(8, 1, 0.5)),
    "Macrophage_24h" = c(rnorm(8, -1, 0.5), rnorm(8, -2, 0.5), rnorm(8, 1.5, 0.5))
  )

  # Create data frame
  expression_matrix <- do.call(cbind, conditions)
  rownames(expression_matrix) <- all_proteins

  # Add p-values (mock)
  p_values <- matrix(runif(length(all_proteins) * (length(conditions) - 1), 0, 0.1),
                     nrow = length(all_proteins))
  colnames(p_values) <- names(conditions)[-1]
  rownames(p_values) <- all_proteins

  # Create mutant comparison data
  mutant_data <- data.frame(
    protein = rep(all_proteins, 4),
    strain = rep(c("WT", "ΔmglA", "ΔqseC", "ΔqseB"), each = length(all_proteins)),
    expression = c(
      rnorm(length(all_proteins), 0, 0.5),  # WT
      c(rnorm(8, -3, 0.5), rnorm(8, -4, 0.5), rnorm(8, -1, 0.5)),  # ΔmglA
      c(rnorm(8, -2, 0.5), rnorm(8, -1, 0.5), rnorm(8, 0, 0.5)),   # ΔqseC
      c(rnorm(8, -2.5, 0.5), rnorm(8, -3, 0.5), rnorm(8, -0.5, 0.5))  # ΔqseB
    ),
    category = rep(c(rep("TCS", 8), rep("FPI", 8), rep("Stress", 8)), 4)
  )

  return(list(
    expression = expression_matrix,
    p_values = p_values,
    mutant_data = mutant_data,
    proteins = all_proteins
  ))
}

# ============================================================================
# FIGURE 1: COMPREHENSIVE TCS REGULATION HEATMAP
# ============================================================================

create_tcs_regulation_heatmap <- function(data) {

  # Filter for TCS proteins only
  tcs_rows <- c("QseB/PmrA", "QseC", "KdpE", "KdpD", "BfpR", "BfpK", "MglA", "SspA")
  tcs_matrix <- data$expression[tcs_rows, ]

  # Create annotation for conditions
  condition_anno <- HeatmapAnnotation(
    Condition = c("Baseline", "Oxidative", "Oxidative", "Acid", "Alkaline",
                  "Heat", "Iron", "Infection", "Infection"),
    Time = c("0h", "2h", "8h", "2h", "2h", "2h", "2h", "2h", "24h"),
    col = list(
      Condition = c("Baseline" = "gray", "Oxidative" = "red", "Acid" = "orange",
                   "Alkaline" = "blue", "Heat" = "purple", "Iron" = "brown",
                   "Infection" = "green"),
      Time = c("0h" = "white", "2h" = "lightblue", "8h" = "blue", "24h" = "darkblue")
    )
  )

  # Create annotation for proteins
  protein_anno <- rowAnnotation(
    Type = c(rep("Response Regulator", 1), rep("Sensor Kinase", 1),
             rep("Response Regulator", 1), rep("Sensor Kinase", 1),
             rep("Response Regulator", 1), rep("Sensor Kinase", 1),
             rep("Global Regulator", 2)),
    Strain = c(rep("All", 2), rep("Environmental", 4), rep("All", 2)),
    col = list(
      Type = c("Response Regulator" = "#E41A1C", "Sensor Kinase" = "#377EB8",
               "Global Regulator" = "#4DAF4A"),
      Strain = c("All" = "black", "Environmental" = "gray50")
    ),
    width = unit(c(0.5, 0.5), "cm")
  )

  # Create the heatmap
  p1 <- Heatmap(
    tcs_matrix,
    name = "Log2 FC",
    col = colorRamp2(c(-4, -2, 0, 2, 4),
                     c("#2166AC", "#4393C3", "white", "#F4A582", "#B2182B")),
    top_annotation = condition_anno,
    left_annotation = protein_anno,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_rot = 45,
    column_title = "TCS Regulation Across Stress Conditions",
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_title = "Two-Component System Proteins",
    row_dend_width = unit(2, "cm"),
    width = ncol(tcs_matrix)*unit(0.8, "cm"),
    height = nrow(tcs_matrix)*unit(0.8, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      # Add significance stars
      if (j > 1) {  # Skip control column
        pval <- data$p_values[rownames(tcs_matrix)[i], colnames(tcs_matrix)[j]]
        if (pval < 0.001) {
          grid.text("***", x, y, gp = gpar(fontsize = 8))
        } else if (pval < 0.01) {
          grid.text("**", x, y, gp = gpar(fontsize = 8))
        } else if (pval < 0.05) {
          grid.text("*", x, y, gp = gpar(fontsize = 8))
        }
      }
    }
  )

  # Draw the heatmap
  draw(p1, heatmap_legend_side = "right")
}

# ============================================================================
# FIGURE 2: WILD-TYPE VS MUTANT STRAIN COMPARISON
# ============================================================================

create_mutant_comparison_plot <- function(data) {

  # Filter for TCS proteins
  tcs_proteins <- c("QseB/PmrA", "QseC", "KdpE", "KdpD", "BfpR", "BfpK", "MglA", "SspA")

  mutant_tcs <- data$mutant_data %>%
    filter(protein %in% tcs_proteins) %>%
    mutate(protein = factor(protein, levels = tcs_proteins))

  # Create grouped barplot
  p2 <- ggplot(mutant_tcs, aes(x = protein, y = expression, fill = strain)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = c("WT" = "#4DAF4A", "ΔmglA" = "#E41A1C",
                                "ΔqseC" = "#377EB8", "ΔqseB" = "#FF7F00")) +
    labs(
      title = "TCS Expression in Wild-Type vs Mutant Strains",
      subtitle = "Baseline expression levels in F. novicida",
      x = "Two-Component System Proteins",
      y = "Log2 Expression (relative to WT)",
      fill = "Strain"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "top"
    ) +
    facet_wrap(~ category, scales = "free_x", ncol = 3) +
    coord_cartesian(ylim = c(-5, 2))

  return(p2)
}

# ============================================================================
# FIGURE 3: TEMPORAL DYNAMICS DURING INFECTION
# ============================================================================

create_temporal_dynamics_plot <- function(data) {

  # Generate time course data
  time_points <- c(0, 0.5, 1, 2, 4, 8, 12, 24)
  proteins_of_interest <- c("QseB/PmrA", "QseC", "MglA", "IglA", "IglB", "IglC")

  # Simulate temporal expression patterns
  temporal_data <- expand.grid(
    time = time_points,
    protein = proteins_of_interest,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      expression = case_when(
        protein %in% c("QseC", "QseB/PmrA") & time <= 2 ~ 2 + time * 1.5 + rnorm(n(), 0, 0.2),
        protein %in% c("QseC", "QseB/PmrA") & time > 2 ~ 5 - (time - 2) * 0.5 + rnorm(n(), 0, 0.2),
        protein == "MglA" ~ 0.5 * time + rnorm(n(), 0, 0.2),
        protein %in% c("IglA", "IglB", "IglC") & time <= 4 ~ 3 + time * 1.5 + rnorm(n(), 0, 0.2),
        protein %in% c("IglA", "IglB", "IglC") & time > 4 ~ 9 - time * 0.8 + rnorm(n(), 0, 0.2),
        TRUE ~ rnorm(n(), 0, 0.5)
      ),
      category = case_when(
        protein %in% c("QseB/PmrA", "QseC") ~ "TCS",
        protein == "MglA" ~ "Global Regulator",
        TRUE ~ "FPI Genes"
      ),
      phase = case_when(
        time <= 4 ~ "Early (Phagosomal)",
        TRUE ~ "Late (Cytosolic)"
      )
    )

  # Create multi-panel temporal plot
  p3 <- ggplot(temporal_data, aes(x = time, y = expression, color = protein)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_vline(xintercept = 4, linetype = "dashed", alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "yellow") +
    annotate("rect", xmin = 4, xmax = 24, ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "blue") +
    annotate("text", x = 2, y = 8, label = "Phagosomal\nEscape",
             size = 3, fontface = "italic") +
    annotate("text", x = 14, y = 8, label = "Cytosolic\nReplication",
             size = 3, fontface = "italic") +
    scale_color_manual(values = c(
      "QseB/PmrA" = "#E41A1C", "QseC" = "#377EB8",
      "MglA" = "#4DAF4A", "IglA" = "#984EA3",
      "IglB" = "#FF7F00", "IglC" = "#A65628"
    )) +
    labs(
      title = "Temporal Dynamics of TCS and Virulence Factors During Infection",
      subtitle = "F. novicida in macrophages",
      x = "Time Post-Infection (hours)",
      y = "Log2 Fold Change (relative to 0h)",
      color = "Protein"
    ) +
    facet_wrap(~ category, ncol = 1, scales = "free_y") +
    theme(
      strip.background = element_rect(fill = "gray90"),
      legend.position = "right"
    )

  return(p3)
}

# ============================================================================
# FIGURE 4: VOLCANO PLOTS FOR KEY CONDITIONS
# ============================================================================

create_volcano_plot_grid <- function(data) {

  # Create volcano data for multiple conditions
  conditions <- c("Oxidative_8h", "pH5.5", "Macrophage_24h", "Iron_limit")

  volcano_plots <- list()

  for (cond_idx in 1:length(conditions)) {
    cond <- conditions[cond_idx]
    # Create volcano data
    volcano_data <- data.frame(
      protein = rownames(data$expression),
      log2FC = data$expression[, cond],
      pvalue = data$p_values[, gsub("_", "_", cond)],
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        padj = p.adjust(pvalue, method = "BH"),
        significance = case_when(
          padj < 0.05 & abs(log2FC) > 1.5 ~ "Significant",
          TRUE ~ "Not Significant"
        ),
        label = ifelse(
          padj < 0.01 & abs(log2FC) > 2 &
          protein %in% c("QseB/PmrA", "QseC", "MglA", "IglA", "KatG"),
          protein, NA
        )
      )

    # Create individual volcano plot
    p <- ggplot(volcano_data, aes(x = log2FC, y = -log10(padj))) +
      geom_point(aes(color = significance), alpha = 0.6, size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", alpha = 0.5) +
      geom_text_repel(aes(label = label), size = 3, max.overlaps = 15) +
      scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray50")) +
      labs(
        title = gsub("_", " ", cond),
        x = "Log2 Fold Change",
        y = "-Log10 (Adjusted P-value)"
      ) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 12)
      ) +
      xlim(c(-6, 6))

    volcano_plots[[cond]] <- p
  }

  # Arrange volcano plots in grid
  p4 <- plot_grid(
    volcano_plots[[1]], volcano_plots[[2]],
    volcano_plots[[3]], volcano_plots[[4]],
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 14
  )

  # Add overall title
  title <- ggdraw() +
    draw_label("TCS Differential Expression Under Key Stress Conditions",
               fontface = 'bold', size = 16)

  p4_final <- plot_grid(title, p4, ncol = 1, rel_heights = c(0.05, 1))

  return(p4_final)
}

# ============================================================================
# FIGURE 5: REGULATORY CASCADE NETWORK
# ============================================================================

create_regulatory_network <- function() {

  # Define network edges
  edges <- data.frame(
    from = c("Membrane Stress", "Membrane Stress", "QseC", "KdpD",
             "QseB/PmrA", "QseB/PmrA", "MglA", "MglA", "MglA",
             "FPI", "FPI", "FPI"),
    to = c("QseC", "KdpD", "QseB/PmrA", "QseB/PmrA",
           "MglA", "FPI", "KatG", "AhpC", "OxyR",
           "IglA", "IglB", "VgrG"),
    interaction = c("disrupts", "disrupts", "phosphorylates", "phosphorylates",
                   "activates", "induces", "induces", "induces", "induces",
                   "encodes", "encodes", "encodes"),
    strength = c(3, 2, 3, 2, 3, 3, 2, 2, 2, 1, 1, 1)
  )

  # Define node properties
  nodes <- data.frame(
    name = unique(c(edges$from, edges$to)),
    type = c("Stress", "Sensor", "Sensor", "Regulator", "Global",
             "Virulence", "Stress Response", "Stress Response", "Stress Response",
             "Effector", "Effector", "Effector"),
    expression = c(0, -2, -1.5, -3, -2.5, -4, 2, 1.5, 1.8, -3, -2.5, -3.5)
  )

  # Create network object
  net <- network(edges[, 1:2], directed = TRUE)

  # Add edge attributes
  set.edge.attribute(net, "interaction", edges$interaction)
  set.edge.attribute(net, "strength", edges$strength)

  # Create ggnetwork data
  net_data <- ggnetwork(net, layout = "kamadakawai")

  # Merge node attributes
  net_data <- merge(net_data, nodes, by.x = "vertex.names", by.y = "name", all.x = TRUE)

  # Create network plot
  p5 <- ggplot(net_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size = strength, alpha = strength),
               arrow = arrow(length = unit(6, "pt"), type = "closed"),
               color = "gray40") +
    geom_nodes(aes(color = expression, shape = type), size = 8) +
    geom_nodelabel_repel(aes(label = vertex.names), size = 3) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "Expression\n(Log2 FC)") +
    scale_shape_manual(values = c("Stress" = 15, "Sensor" = 17, "Regulator" = 16,
                                 "Global" = 18, "Virulence" = 19,
                                 "Stress Response" = 8, "Effector" = 13)) +
    labs(
      title = "TCS Regulatory Cascade Under Membrane Stress",
      subtitle = "Orphan system vulnerability in F. novicida",
      shape = "Component Type",
      size = "Interaction\nStrength",
      alpha = "Interaction\nStrength"
    ) +
    theme_blank() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    )

  return(p5)
}

# ============================================================================
# FIGURE 6: COMPARATIVE ANALYSIS ACROSS SPECIES/STRAINS
# ============================================================================

create_species_comparison <- function() {

  # Create comparison data
  comparison_data <- expand.grid(
    protein = c("QseB/PmrA", "QseC", "KdpE", "KdpD", "BfpR", "BfpK"),
    strain = c("F. novicida", "F. tularensis SchuS4", "F. tularensis LVS"),
    condition = c("Control", "Oxidative", "pH Stress", "Iron Limit")
  ) %>%
    mutate(
      functional = case_when(
        strain == "F. novicida" ~ TRUE,
        strain == "F. tularensis SchuS4" & protein %in% c("KdpE", "BfpK") ~ FALSE,
        strain == "F. tularensis LVS" & protein %in% c("KdpE", "KdpD", "BfpR", "BfpK") ~ FALSE,
        TRUE ~ TRUE
      ),
      expression = ifelse(
        functional,
        rnorm(n(), mean = ifelse(condition == "Control", 0,
                                ifelse(condition == "Oxidative", 2,
                                     ifelse(condition == "pH Stress", -1, 1.5))),
              sd = 0.5),
        NA
      )
    )

  # Create heatmap-style comparison
  p6 <- ggplot(comparison_data, aes(x = protein, y = strain, fill = expression)) +
    geom_tile(aes(alpha = functional), color = "white", size = 0.5) +
    geom_text(aes(label = ifelse(is.na(expression), "Ψ", round(expression, 1))),
              size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, na.value = "gray90",
                        name = "Log2 FC") +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1),
                      guide = "none") +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      title = "TCS Functional Capacity Across Francisella Strains",
      subtitle = "Ψ indicates pseudogene or absent",
      x = "Two-Component System Gene",
      y = "Strain"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "gray90"),
      panel.background = element_rect(fill = "gray95")
    )

  return(p6)
}

# ============================================================================
# FIGURE 7: SUMMARY STATISTICS AND KEY FINDINGS
# ============================================================================

create_summary_statistics_plot <- function(data) {

  # Calculate summary statistics
  summary_stats <- data.frame(
    condition = colnames(data$expression)[-1],
    n_upregulated = colSums(data$expression[, -1] > 1.5),
    n_downregulated = colSums(data$expression[, -1] < -1.5),
    mean_fc_tcs = colMeans(data$expression[1:8, -1]),
    mean_fc_fpi = colMeans(data$expression[9:16, -1])
  ) %>%
    pivot_longer(cols = -condition, names_to = "metric", values_to = "value") %>%
    mutate(
      category = case_when(
        grepl("upregulated", metric) ~ "Upregulated Proteins",
        grepl("downregulated", metric) ~ "Downregulated Proteins",
        grepl("tcs", metric) ~ "Mean TCS FC",
        grepl("fpi", metric) ~ "Mean FPI FC"
      )
    )

  # Create summary plots
  p7a <- ggplot(filter(summary_stats, grepl("regulated", metric)),
               aes(x = condition, y = value, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Upregulated Proteins" = "#E41A1C",
                                "Downregulated Proteins" = "#377EB8")) +
    labs(
      title = "Number of Regulated Proteins by Condition",
      x = "Condition",
      y = "Number of Proteins",
      fill = "Direction"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p7b <- ggplot(filter(summary_stats, grepl("Mean", category)),
               aes(x = condition, y = value, color = category, group = category)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("Mean TCS FC" = "#4DAF4A",
                                 "Mean FPI FC" = "#FF7F00")) +
    labs(
      title = "Average Expression Changes: TCS vs FPI Genes",
      x = "Condition",
      y = "Mean Log2 Fold Change",
      color = "Gene Category"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Combine plots
  p7 <- plot_grid(p7a, p7b, ncol = 1, labels = c("A", "B"))

  return(p7)
}

# ============================================================================
# MAIN EXECUTION: GENERATE ALL FIGURES
# ============================================================================

# Generate mock data (replace with actual ProteomeXchange data loading)
tcs_data <- generate_mock_tcs_data()

# Create output directory
dir.create("figures", showWarnings = FALSE)

# Generate Figure 1: Comprehensive TCS Heatmap
pdf("figures/Fig1_TCS_Regulation_Heatmap.pdf", width = 12, height = 8)
create_tcs_regulation_heatmap(tcs_data)
dev.off()

# Generate Figure 2: Wild-type vs Mutant Comparison
pdf("figures/Fig2_Mutant_Comparison.pdf", width = 14, height = 6)
print(create_mutant_comparison_plot(tcs_data))
dev.off()

# Generate Figure 3: Temporal Dynamics
pdf("figures/Fig3_Temporal_Dynamics.pdf", width = 10, height = 10)
print(create_temporal_dynamics_plot(tcs_data))
dev.off()

# Generate Figure 4: Volcano Plot Grid
pdf("figures/Fig4_Volcano_Plots.pdf", width = 12, height = 10)
print(create_volcano_plot_grid(tcs_data))
dev.off()

# Generate Figure 5: Regulatory Network
pdf("figures/Fig5_Regulatory_Network.pdf", width = 10, height = 8)
print(create_regulatory_network())
dev.off()

# Generate Figure 6: Species Comparison
pdf("figures/Fig6_Species_Comparison.pdf", width = 12, height = 8)
print(create_species_comparison())
dev.off()

# Generate Figure 7: Summary Statistics
pdf("figures/Fig7_Summary_Statistics.pdf", width = 10, height = 10)
print(create_summary_statistics_plot(tcs_data))
dev.off()

# ============================================================================
# COMPOSITE FIGURE FOR PUBLICATION
# ============================================================================

create_composite_figure <- function() {

  # Load individual plots
  p1 <- create_mutant_comparison_plot(tcs_data)
  p2 <- create_temporal_dynamics_plot(tcs_data)
  p3 <- create_regulatory_network()
  p4 <- create_species_comparison()

  # Create composite
  composite <- plot_grid(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 16,
    rel_widths = c(1, 1),
    rel_heights = c(1, 1)
  )

  # Add title
  title <- ggdraw() +
    draw_label("Two-Component System Regulation in Francisella novicida",
               fontface = 'bold', size = 18)

  final_composite <- plot_grid(title, composite, ncol = 1,
                              rel_heights = c(0.05, 1))

  return(final_composite)
}

# Generate composite figure
pdf("figures/Fig_Composite_TCS_Analysis.pdf", width = 16, height = 14)
print(create_composite_figure())
dev.off()

# ============================================================================
# SUPPLEMENTARY FIGURES
# ============================================================================

# Generate supplementary table of all significant TCS changes
create_supplementary_table <- function(data) {

  # Extract significant changes
  sig_changes <- data.frame()

  for (col in colnames(data$expression)[-1]) {
    temp <- data.frame(
      protein = rownames(data$expression),
      condition = col,
      log2FC = data$expression[, col],
      pvalue = data$p_values[, col],
      padj = p.adjust(data$p_values[, col], method = "BH")
    ) %>%
      filter(padj < 0.05) %>%
      arrange(padj)

    sig_changes <- rbind(sig_changes, temp)
  }

  # Filter for TCS proteins
  tcs_sig <- sig_changes %>%
    filter(grepl("Qse|Kdp|Bfp|Mgl|Ssp", protein)) %>%
    mutate(
      direction = ifelse(log2FC > 0, "UP", "DOWN"),
      magnitude = case_when(
        abs(log2FC) > 3 ~ "Strong",
        abs(log2FC) > 1.5 ~ "Moderate",
        TRUE ~ "Mild"
      )
    )

  # Save to CSV
  write.csv(tcs_sig, "figures/Table_S1_TCS_Significant_Changes.csv", row.names = FALSE)

  return(tcs_sig)
}

# Generate supplementary table
supp_table <- create_supplementary_table(tcs_data)

message("All figures generated successfully!")
message("Check the 'figures' directory for output files")

# ============================================================================
# FIGURE LEGENDS FOR MANUSCRIPT
# ============================================================================

figure_legends <- list(
  Fig1 = "Figure 1. Comprehensive regulation of two-component system proteins across multiple stress conditions in F. novicida. Heatmap showing log2 fold changes relative to control conditions. Hierarchical clustering reveals coordinated regulation patterns. Asterisks indicate statistical significance (*p<0.05, **p<0.01, ***p<0.001).",

  Fig2 = "Figure 2. Comparative expression of TCS proteins in wild-type versus isogenic mutant strains. Deletion of global regulators (ΔmglA) shows more severe effects than individual TCS deletions (ΔqseC, ΔqseB), supporting the orphan system vulnerability hypothesis.",

  Fig3 = "Figure 3. Temporal dynamics of TCS and virulence factor expression during macrophage infection. Early upregulation of TCS proteins (0-4h) corresponds to phagosomal escape phase, while FPI gene expression peaks during early cytosolic replication (2-8h) then declines.",

  Fig4 = "Figure 4. Volcano plots showing differential expression under key stress conditions. TCS proteins show condition-specific regulation patterns, with oxidative stress and macrophage infection causing the most dramatic changes.",

  Fig5 = "Figure 5. Regulatory network illustrating the cascade effects of membrane stress on the orphan TCS system. Network analysis reveals the vulnerability of having limited sensor kinases (QseC, KdpD) controlling critical virulence regulators.",

  Fig6 = "Figure 6. Strain-specific TCS functional capacity across Francisella species. Environmental strains (F. novicida) maintain complete TCS arsenals while pathogenic strains show progressive gene loss, creating differential susceptibility to stress.",

  Fig7 = "Figure 7. Summary statistics of proteomic changes. (A) Number of significantly regulated proteins per condition. (B) Average expression changes comparing TCS versus FPI genes across conditions, revealing inverse regulation patterns."
)

# Save figure legends
writeLines(unlist(figure_legends), "figures/Figure_Legends.txt")
