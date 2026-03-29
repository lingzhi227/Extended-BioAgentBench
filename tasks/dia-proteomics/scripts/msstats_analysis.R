#!/usr/bin/env Rscript
# MSstats analysis for DIA proteomics data
# Input: TRIC-aligned feature file + sample sheet
# Output: differential expression results

suppressPackageStartupMessages({
  library(MSstats)
})

args <- commandArgs(trailingOnly = TRUE)
aligned_file <- args[1]
sample_sheet_file <- args[2]
output_dir <- args[3]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== MSstats Analysis ===\n")

# Read aligned features from TRIC
aligned <- read.delim(aligned_file, stringsAsFactors = FALSE)
cat("Aligned features:", nrow(aligned), "\n")

# Read sample sheet
sample_sheet <- read.delim(sample_sheet_file, stringsAsFactors = FALSE)
cat("Sample sheet:\n")
print(sample_sheet)

# Prepare MSstats input format
# The aligned file from TRIC has columns like:
# transition_group_id, run_id, Intensity, etc.

# Map run names to conditions from sample sheet
tryCatch({
  # Check available columns
  cat("Aligned columns:", paste(head(names(aligned), 20), collapse=", "), "\n")

  # Try to build MSstats-compatible input
  # OpenSWATH/TRIC output needs conversion
  if ("ProteinName" %in% names(aligned) && "filename" %in% names(aligned)) {
    # Create condition mapping
    run_to_condition <- setNames(sample_sheet$Condition, sample_sheet$Sample)

    # Build MSstats input
    msstats_input <- data.frame(
      ProteinName = aligned$ProteinName,
      PeptideSequence = if ("FullPeptideName" %in% names(aligned)) aligned$FullPeptideName else aligned$Sequence,
      PrecursorCharge = if ("Charge" %in% names(aligned)) aligned$Charge else 2,
      FragmentIon = if ("aggr_Fragment_Annotation" %in% names(aligned)) aligned$aggr_Fragment_Annotation else "y",
      ProductCharge = 1,
      IsotopeLabelType = "L",
      Condition = sapply(aligned$filename, function(fn) {
        # Match filename to sample sheet
        matched <- which(sapply(sample_sheet$Sample, function(s) grepl(s, fn, fixed=TRUE)))
        if (length(matched) > 0) sample_sheet$Condition[matched[1]] else "Unknown"
      }),
      BioReplicate = aligned$filename,
      Run = aligned$filename,
      Intensity = if ("Intensity" %in% names(aligned)) aligned$Intensity else
                  if ("m_score" %in% names(aligned)) 10^6 else 0,
      stringsAsFactors = FALSE
    )

    cat("MSstats input rows:", nrow(msstats_input), "\n")
    cat("Conditions:", paste(unique(msstats_input$Condition), collapse=", "), "\n")

    # Run MSstats
    processed <- dataProcess(msstats_input, logTrans = 2, normalization = "equalizeMedians")

    # Check if we have 2+ conditions
    conditions <- unique(msstats_input$Condition)
    conditions <- conditions[conditions != "Unknown"]

    if (length(conditions) >= 2) {
      contrast_matrix <- matrix(c(1, -1), nrow=1)
      colnames(contrast_matrix) <- conditions[1:2]
      rownames(contrast_matrix) <- paste(conditions[1], "vs", conditions[2])

      result <- groupComparison(contrast.matrix = contrast_matrix, data = processed)
      comparison <- result$ComparisonResult

      write.csv(comparison, file.path(output_dir, "msstats_results.csv"), row.names = FALSE)
      cat("DE proteins (adj.pvalue < 0.05):", sum(comparison$adj.pvalue < 0.05, na.rm=TRUE), "\n")
    } else {
      cat("Only one condition found, skipping differential analysis\n")
      # Write protein-level summary instead
      protein_summary <- data.frame(
        Protein = unique(msstats_input$ProteinName),
        adj.pvalue = NA
      )
      write.csv(protein_summary, file.path(output_dir, "msstats_results.csv"), row.names = FALSE)
    }
  } else {
    cat("Required columns not found in aligned file\n")
    cat("Available columns:", paste(names(aligned), collapse=", "), "\n")
    # Write placeholder
    write.csv(data.frame(Protein=character(0), adj.pvalue=numeric(0)),
              file.path(output_dir, "msstats_results.csv"), row.names = FALSE)
  }
}, error = function(e) {
  cat("MSstats error:", conditionMessage(e), "\n")
  write.csv(data.frame(Protein=character(0), adj.pvalue=numeric(0)),
            file.path(output_dir, "msstats_results.csv"), row.names = FALSE)
})

cat("=== MSstats complete ===\n")
