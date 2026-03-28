#!/usr/bin/env Rscript
# XCMS processing: read mzML -> peak detection -> grouping -> alignment -> fill
# Outputs: peak_matrix.csv, feature_defs.csv, xcms_summary.csv, xdata.rds

suppressPackageStartupMessages({
  library(MSnbase)
  library(xcms)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_dir <- args[2]
threads <- as.integer(args[3])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Register parallel backend
if (threads > 1) {
  register(bpstart(MulticoreParam(threads)))
}

cat("=== Step 1: Reading mzML files ===\n")
mzml_files <- sort(list.files(data_dir, pattern = "\\.mzML$", full.names = TRUE))
cat("Found", length(mzml_files), "mzML files\n")
raw_data <- readMSData(mzml_files, mode = "onDisk")

cat("=== Step 2: Peak detection (MatchedFilter for GC-MS) ===\n")
mfp <- MatchedFilterParam(fwhm = 10, snthresh = 5, binSize = 0.5,
                           steps = 2, mzdiff = 0.5)
xdata <- findChromPeaks(raw_data, param = mfp)
n_peaks_raw <- nrow(chromPeaks(xdata))
cat("Raw peaks detected:", n_peaks_raw, "\n")

cat("=== Step 3: Initial peak grouping ===\n")
sample_groups <- rep("sample", length(mzml_files))
pdp <- PeakDensityParam(sampleGroups = sample_groups,
                         bw = 5, minFraction = 0.5,
                         binSize = 0.25)
xdata <- groupChromPeaks(xdata, param = pdp)
n_features_pre <- nrow(featureDefinitions(xdata))
cat("Features before alignment:", n_features_pre, "\n")

cat("=== Step 4: Retention time alignment ===\n")
pgp <- PeakGroupsParam(minFraction = 0.8, smooth = "loess",
                         span = 0.4, family = "gaussian")
xdata <- adjustRtime(xdata, param = pgp)

cat("=== Step 5: Re-grouping after alignment ===\n")
xdata <- groupChromPeaks(xdata, param = pdp)
n_features_post <- nrow(featureDefinitions(xdata))
cat("Features after alignment:", n_features_post, "\n")

cat("=== Step 6: Fill missing peaks ===\n")
xdata <- fillChromPeaks(xdata)

# Export peak matrix (samples x features, intensity values)
feat_vals <- featureValues(xdata, value = "into")
feat_defs <- featureDefinitions(xdata)

# Convert list columns to character for CSV export
feat_defs_df <- as.data.frame(feat_defs)
for (col in names(feat_defs_df)) {
  if (is.list(feat_defs_df[[col]])) {
    feat_defs_df[[col]] <- sapply(feat_defs_df[[col]], function(x) paste(x, collapse = ";"))
  }
}

write.csv(feat_vals, file.path(output_dir, "peak_matrix.csv"))
write.csv(feat_defs_df, file.path(output_dir, "feature_defs.csv"))

# Save R object for downstream use
saveRDS(xdata, file.path(output_dir, "xdata.rds"))

# Summary stats
n_features <- nrow(feat_defs)
n_samples <- ncol(feat_vals)
rt_range <- paste0(round(min(feat_defs$rtmin, na.rm=TRUE), 2), "-",
                    round(max(feat_defs$rtmax, na.rm=TRUE), 2))
mz_range <- paste0(round(min(feat_defs$mzmin, na.rm=TRUE), 2), "-",
                    round(max(feat_defs$mzmax, na.rm=TRUE), 2))

# Feature intensity stats
mean_intensities <- rowMeans(feat_vals, na.rm = TRUE)
fill_rate <- sum(!is.na(feat_vals)) / length(feat_vals) * 100

summary_df <- data.frame(
  metric = c("total_peaks_detected", "features_before_alignment",
             "features_after_alignment", "samples_processed",
             "rt_range_minutes", "mz_range", "fill_rate_percent"),
  value = c(n_peaks_raw, n_features_pre, n_features_post,
            n_samples, rt_range, mz_range, round(fill_rate, 1))
)
write.csv(summary_df, file.path(output_dir, "xcms_summary.csv"), row.names = FALSE)
cat("=== XCMS processing complete ===\n")
cat("Features:", n_features, "| Samples:", n_samples, "\n")
