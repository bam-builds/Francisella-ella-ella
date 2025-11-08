#!/usr/bin/env Rscript
# ============================================================
# MaxQuant Data Preparation for DEqMS Analysis
# ============================================================
# Converts MaxQuant proteinGroups.txt and evidence.txt files
# into format suitable for DEqMS analysis
#
# Author: DEqMS TCS Analysis Pipeline
# Date: 2025
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

cat("=============================================================\n")
cat("MaxQuant Data Preparation for DEqMS\n")
cat("=============================================================\n\n")

# ============================================================
# CONFIGURATION
# ============================================================

# Input files (MaxQuant output)
MAXQUANT_DIR <- "data/raw/maxquant"
PROTEIN_GROUPS_FILE <- file.path(MAXQUANT_DIR, "proteinGroups.txt")
EVIDENCE_FILE <- file.path(MAXQUANT_DIR, "evidence.txt")

# Output files
OUTPUT_DIR <- "data/processed"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_PROTEIN_FILE <- file.path(OUTPUT_DIR, "PXD035145_proteinGroups.csv")
OUTPUT_PSM_FILE <- file.path(OUTPUT_DIR, "PXD035145_PSM_counts.csv")

# ============================================================
# FUNCTION: Extract Intensity Columns
# ============================================================

#' Extract LFQ or intensity columns from MaxQuant output
#'
#' @param df Data frame from proteinGroups.txt
#' @param intensity_type Type of intensity ("LFQ" or "Intensity")
#' @return Matrix of intensity values
extract_intensities <- function(df, intensity_type = "LFQ") {

  if (intensity_type == "LFQ") {
    intensity_cols <- grep("^LFQ\\.intensity\\.", colnames(df), value = TRUE)
  } else {
    intensity_cols <- grep("^Intensity\\.", colnames(df), value = TRUE)
    # Exclude LFQ columns if present
    intensity_cols <- intensity_cols[!grepl("LFQ", intensity_cols)]
  }

  if (length(intensity_cols) == 0) {
    stop("No ", intensity_type, " intensity columns found in data")
  }

  cat("  Found", length(intensity_cols), intensity_type, "intensity columns\n")

  # Extract intensity matrix
  intensity_mat <- as.matrix(df[, intensity_cols])

  # Clean sample names
  sample_names <- gsub(paste0("^", intensity_type, "\\.intensity\\."), "",
                      colnames(intensity_mat))
  sample_names <- gsub("^Intensity\\.", "", sample_names)
  colnames(intensity_mat) <- sample_names

  # Use gene names as row names (or protein IDs if gene names missing)
  if ("Gene.names" %in% colnames(df)) {
    rownames(intensity_mat) <- ifelse(
      df$Gene.names != "",
      df$Gene.names,
      df$Protein.IDs
    )
  } else if ("Protein.IDs" %in% colnames(df)) {
    rownames(intensity_mat) <- df$Protein.IDs
  } else {
    rownames(intensity_mat) <- paste0("Protein_", 1:nrow(df))
  }

  # Convert zeros to NA (common in MaxQuant output)
  intensity_mat[intensity_mat == 0] <- NA

  return(intensity_mat)
}

# ============================================================
# FUNCTION: Calculate PSM Counts from Evidence
# ============================================================

#' Calculate PSM counts per protein from evidence.txt
#'
#' @param evidence_file Path to evidence.txt
#' @param protein_ids Vector of protein IDs to match
#' @return Data frame with PSM counts per protein
calculate_psm_counts <- function(evidence_file, protein_ids) {

  cat("\nCalculating PSM counts from evidence file...\n")

  if (!file.exists(evidence_file)) {
    warning("Evidence file not found: ", evidence_file)
    warning("Using peptide counts as proxy for PSM counts")
    return(NULL)
  }

  # Read evidence file (can be large, so read selectively)
  evidence <- read.delim(evidence_file, stringsAsFactors = FALSE)

  cat("  Loaded", nrow(evidence), "rows from evidence.txt\n")

  # Count PSMs per protein
  # Note: In MaxQuant, each row in evidence.txt is a peptide-spectrum match
  psm_counts <- evidence %>%
    filter(!is.na(Proteins), Proteins != "") %>%
    group_by(Proteins) %>%
    summarise(PSM_count = n(), .groups = "drop")

  cat("  Calculated PSM counts for", nrow(psm_counts), "proteins\n")

  return(psm_counts)
}

# ============================================================
# MAIN PROCESSING
# ============================================================

cat("Step 1: Loading proteinGroups.txt\n")
cat("-------------------------------------------------------------\n")

if (!file.exists(PROTEIN_GROUPS_FILE)) {
  stop("ERROR: proteinGroups.txt not found at: ", PROTEIN_GROUPS_FILE, "\n",
       "Please ensure MaxQuant output is in: ", MAXQUANT_DIR)
}

# Read protein groups
protein_groups <- read.delim(PROTEIN_GROUPS_FILE, stringsAsFactors = FALSE)
cat("  Loaded", nrow(protein_groups), "protein groups\n\n")

# ============================================================
# Step 2: Filter Contaminants and Reverse Hits
# ============================================================

cat("Step 2: Filtering Contaminants and Decoys\n")
cat("-------------------------------------------------------------\n")

original_count <- nrow(protein_groups)

# Filter reverse hits and contaminants
protein_groups_clean <- protein_groups %>%
  filter(
    # Remove reverse hits (decoys)
    is.na(Reverse) | Reverse != "+",
    # Remove contaminants
    is.na(Potential.contaminant) | Potential.contaminant != "+",
    # Remove proteins only identified by site
    is.na(Only.identified.by.site) | Only.identified.by.site != "+"
  )

cat("  Original proteins:      ", original_count, "\n")
cat("  After filtering:        ", nrow(protein_groups_clean), "\n")
cat("  Removed:                ", original_count - nrow(protein_groups_clean), "\n\n")

# ============================================================
# Step 3: Extract Intensity Matrix
# ============================================================

cat("Step 3: Extracting Intensity Values\n")
cat("-------------------------------------------------------------\n")

# Try LFQ intensities first, fall back to regular intensities
intensity_matrix <- tryCatch({
  extract_intensities(protein_groups_clean, intensity_type = "LFQ")
}, error = function(e) {
  cat("  LFQ intensities not found, using regular intensities\n")
  extract_intensities(protein_groups_clean, intensity_type = "Intensity")
})

cat("  Intensity matrix: ", nrow(intensity_matrix), " proteins × ",
    ncol(intensity_matrix), " samples\n\n")

# ============================================================
# Step 4: Extract or Calculate PSM Counts
# ============================================================

cat("Step 4: PSM Count Extraction\n")
cat("-------------------------------------------------------------\n")

# First, try to get PSM counts from evidence.txt
psm_counts_evidence <- calculate_psm_counts(
  EVIDENCE_FILE,
  protein_groups_clean$Protein.IDs
)

# Create PSM count matrix
if (!is.null(psm_counts_evidence)) {
  # Match PSM counts to protein matrix
  psm_matrix <- data.frame(
    Protein = rownames(intensity_matrix),
    PSM_count = 1,  # Default value
    stringsAsFactors = FALSE
  )

  # Map evidence PSM counts to gene names/protein IDs
  for (i in 1:nrow(psm_matrix)) {
    protein_id <- psm_matrix$Protein[i]

    # Try to find matching entry in evidence PSM counts
    # This requires mapping gene names back to protein IDs
    matching_idx <- which(protein_groups_clean$Gene.names == protein_id |
                         protein_groups_clean$Protein.IDs == protein_id)

    if (length(matching_idx) > 0) {
      original_id <- protein_groups_clean$Protein.IDs[matching_idx[1]]
      psm_match <- psm_counts_evidence %>%
        filter(str_detect(Proteins, fixed(original_id)))

      if (nrow(psm_match) > 0) {
        psm_matrix$PSM_count[i] <- psm_match$PSM_count[1]
      }
    }
  }

} else {
  # Use peptide counts as fallback
  cat("  Using peptide counts as proxy for PSM counts\n")

  peptide_counts <- protein_groups_clean %>%
    select(Protein.IDs, Gene.names, Peptides = Peptides) %>%
    mutate(
      Protein = ifelse(Gene.names != "", Gene.names, Protein.IDs),
      PSM_count = pmax(Peptides, 1)  # Ensure minimum of 1
    ) %>%
    select(Protein, PSM_count)

  # Match to intensity matrix
  psm_matrix <- data.frame(
    Protein = rownames(intensity_matrix),
    PSM_count = 1,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(psm_matrix)) {
    match_idx <- which(peptide_counts$Protein == psm_matrix$Protein[i])
    if (length(match_idx) > 0) {
      psm_matrix$PSM_count[i] <- peptide_counts$PSM_count[match_idx[1]]
    }
  }
}

rownames(psm_matrix) <- psm_matrix$Protein
psm_matrix <- psm_matrix[, "PSM_count", drop = FALSE]

cat("  PSM count summary:\n")
print(summary(psm_matrix$PSM_count))
cat("\n")

# ============================================================
# Step 5: Save Processed Data
# ============================================================

cat("Step 5: Saving Processed Data\n")
cat("-------------------------------------------------------------\n")

# Save intensity matrix
write.csv(intensity_matrix, OUTPUT_PROTEIN_FILE, row.names = TRUE)
cat("  Saved protein intensities: ", OUTPUT_PROTEIN_FILE, "\n")

# Save PSM counts
write.csv(psm_matrix, OUTPUT_PSM_FILE, row.names = TRUE)
cat("  Saved PSM counts: ", OUTPUT_PSM_FILE, "\n\n")

# ============================================================
# Step 6: Summary Report
# ============================================================

cat("=============================================================\n")
cat("DATA PREPARATION SUMMARY\n")
cat("=============================================================\n\n")

cat("Input:\n")
cat("  - proteinGroups.txt:    ", original_count, " entries\n")
cat("  - evidence.txt:         ", ifelse(!is.null(psm_counts_evidence), "processed", "not found"), "\n\n")

cat("Output:\n")
cat("  - Protein matrix:       ", nrow(intensity_matrix), " × ", ncol(intensity_matrix), "\n")
cat("  - PSM counts:           ", nrow(psm_matrix), " proteins\n")
cat("  - Missing values:       ", sprintf("%.1f%%",
    sum(is.na(intensity_matrix)) / length(intensity_matrix) * 100), "\n\n")

cat("Sample names:\n")
cat(" ", paste(colnames(intensity_matrix), collapse = "\n  "), "\n\n")

cat("Ready for DEqMS analysis!\n")
cat("Next step: Run deqms_tcs_analysis.R\n\n")

cat("=============================================================\n")
cat("PREPARATION COMPLETE\n")
cat("=============================================================\n")
