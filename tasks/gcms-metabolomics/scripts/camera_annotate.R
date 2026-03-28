#!/usr/bin/env Rscript
# CAMERA annotation: isotope patterns + adduct groups + pseudospectra
# Input: xdata.rds from xcms_process.R
# Output: camera_annotations.csv, pseudospectra.csv

suppressPackageStartupMessages({
  library(xcms)
  library(CAMERA)
})

args <- commandArgs(trailingOnly = TRUE)
xcms_dir <- args[1]    # directory with xdata.rds
output_dir <- args[2]  # output directory

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Loading XCMS result ===\n")
xdata <- readRDS(file.path(xcms_dir, "xdata.rds"))

# Convert to old-style xcmsSet for CAMERA compatibility
xset <- as(xdata, "xcmsSet")
sampclass(xset) <- rep("sample", length(sampnames(xset)))

cat("=== Running CAMERA annotation ===\n")
# Create annotated object
an <- xsAnnotate(xset)

# Group by RT correlation
an <- groupFWHM(an, perfwhm = 0.6)

# Find isotope patterns
an <- findIsotopes(an, mzabs = 0.5, ppm = 10)

# Find adduct patterns
an <- findAdducts(an, polarity = "positive", ppm = 10)

# Get peak table with annotations
peaklist <- getPeaklist(an)

# Extract annotation summary
n_isotope_groups <- sum(grepl("\\[M\\]", peaklist$isotopes) |
                         grepl("\\[M\\+1\\]", peaklist$isotopes))
n_adduct_annotations <- sum(peaklist$adduct != "")
n_pseudospectra <- max(as.numeric(peaklist$pcgroup), na.rm = TRUE)

cat("Isotope annotations:", n_isotope_groups, "\n")
cat("Adduct annotations:", n_adduct_annotations, "\n")
cat("Pseudospectra groups:", n_pseudospectra, "\n")

# Write full annotated peak list
write.csv(peaklist, file.path(output_dir, "camera_peaklist.csv"), row.names = FALSE)

# Write pseudospectra info (group membership + mz + rt)
ps_info <- data.frame(
  feature_idx = seq_len(nrow(peaklist)),
  pcgroup = peaklist$pcgroup,
  mz = peaklist$mz,
  rt = peaklist$rt,
  isotopes = peaklist$isotopes,
  adduct = peaklist$adduct
)
write.csv(ps_info, file.path(output_dir, "pseudospectra_info.csv"), row.names = FALSE)

# Summary
summary_df <- data.frame(
  metric = c("pseudospectra_count", "isotope_annotations", "adduct_annotations",
             "total_annotated_features"),
  value = c(n_pseudospectra, n_isotope_groups, n_adduct_annotations,
            sum(peaklist$isotopes != "" | peaklist$adduct != ""))
)
write.csv(summary_df, file.path(output_dir, "camera_summary.csv"), row.names = FALSE)
cat("=== CAMERA annotation complete ===\n")
