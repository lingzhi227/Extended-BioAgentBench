#!/usr/bin/env bash
set -euo pipefail

# LC-MS Untargeted Metabolomics Pipeline
# Data: Sacurine dataset (human urine, LTQ-Orbitrap, negative mode)
# Samples: 6 study, 3 QC pools, 3 blanks
#
# DAG Structure (depth=10, convergence=4):
#
#   [study mzML]────[QC mzML]────[blank mzML]
#       │               │             │
#   read mzML       read mzML    read mzML           (Step 1)
#       │               │             │
#       └───────────────┼─────────────┘
#                       │
#                  peak detection                      (Step 2: CONVERGE #1)
#                       │
#              ╱────────┴────────╲
#        RT alignment      peak grouping               (Step 3: parallel)
#              ╲────────┬────────╱
#                       │
#                 fill missing peaks                    (Step 4: CONVERGE #2)
#                       │
#            ╱──────────┼──────────╲
#      isotope/adduct  blank       QC-based             (Step 5: parallel)
#      annotation     subtraction  filtering
#            ╲──────────┼──────────╱
#                       │
#                filtered feature table                 (Step 6: CONVERGE #3)
#                       │
#                ╱──────┴──────╲
#          spectral          statistical                (Step 7: parallel)
#          matching          analysis
#                ╲──────┬──────╱
#                       │
#                 pathway enrichment                    (Step 8: CONVERGE #4)
#                       │
#                   report.csv                          (Step 9-10)

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="${WORKDIR}/data"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}" "${RESULTS}"

# ============================================================================
# Steps 1-8: Run entire R pipeline via XCMS + CAMERA
# ============================================================================
echo "[Pipeline] Running XCMS/CAMERA metabolomics pipeline..."

if [ ! -f "${RESULTS}/report.csv" ]; then

# Write R script to temp file (avoids heredoc escaping issues)
cat > "${OUT}/_pipeline.R" << 'REOF'
library(xcms)
library(MSnbase)

cat("[Step 1] Reading mzML files...\n")
args <- commandArgs(trailingOnly=TRUE)
data_dir <- args[1]
out_dir <- args[2]
results_dir <- args[3]

# Read sample metadata
meta <- read.delim(file.path(data_dir, "sampleMetadata.tsv"), stringsAsFactors=FALSE)

# Identify file paths matching metadata
mzml_files <- file.path(data_dir, paste0(meta$sample_name, ".mzML"))
stopifnot(all(file.exists(mzml_files)))

# Read raw data
raw_data <- readMSData(mzml_files, mode="onDisk", msLevel.=1)

# Assign sample classes
pData(raw_data)$sample_name <- meta$sample_name
pData(raw_data)$sample_type <- meta$sampleType
pData(raw_data)$class <- meta$class
pData(raw_data)$gender <- meta$gender
pData(raw_data)$age <- meta$age
pData(raw_data)$bmi <- meta$bmi
pData(raw_data)$injection_order <- meta$injectionOrder

cat("  Loaded", length(mzml_files), "files\n")
cat("  Sample types:", paste(table(meta$sampleType), collapse=", "), "\n")

# ============================================================================
# Step 2: Peak detection with CentWave — CONVERGE #1 (all sample types)
# ============================================================================
cat("[Step 2] Peak detection (CentWave)...\n")
cwp <- CentWaveParam(
    ppm = 5,
    peakwidth = c(5, 20),
    snthresh = 5,
    mzdiff = 0.01,
    prefilter = c(3, 500),
    noise = 500,
    integrate = 1
)
xdata <- findChromPeaks(raw_data, param=cwp)
cat("  Peaks found:", nrow(chromPeaks(xdata)), "\n")

# ============================================================================
# Step 3a: Retention time alignment (obiwarp)
# ============================================================================
cat("[Step 3a] RT alignment (obiwarp)...\n")
obi <- ObiwarpParam(binSize=0.5)
xdata <- adjustRtime(xdata, param=obi)
cat("  RT correction applied\n")

# ============================================================================
# Step 3b: Peak grouping (peak density)
# ============================================================================
cat("[Step 3b] Peak grouping (density)...\n")
pdp <- PeakDensityParam(
    sampleGroups = pData(xdata)$sample_type,
    minFraction = 0.5,
    bw = 5,
    binSize = 0.01
)
xdata <- groupChromPeaks(xdata, param=pdp)
n_features_grouped <- nrow(featureDefinitions(xdata))
cat("  Features after grouping:", n_features_grouped, "\n")

# ============================================================================
# Step 4: Fill missing peaks — CONVERGE #2 (aligned + grouped)
# ============================================================================
cat("[Step 4] Filling missing peaks (convergence #2)...\n")
xdata <- fillChromPeaks(xdata, param=ChromPeakAreaParam())
n_features_filled <- nrow(featureDefinitions(xdata))
cat("  Features after fill:", n_features_filled, "\n")

# Extract feature values
fvals <- featureValues(xdata, value="into", method="maxint")
fdefs <- featureDefinitions(xdata)

# ============================================================================
# Step 5a: Isotope/adduct annotation (CAMERA)
# ============================================================================
cat("[Step 5a] Isotope and adduct annotation (CAMERA)...\n")
suppressMessages(library(CAMERA))
xset <- as(xdata, "xcmsSet")
an <- xsAnnotate(xset)
an <- groupFWHM(an, perfwhm=0.6)
an <- findIsotopes(an, mzabs=0.01, ppm=5)
an <- findAdducts(an, polarity="negative")
peaklist <- getPeaklist(an)

n_isotopes <- sum(peaklist$isotopes != "", na.rm=TRUE)
n_adducts <- sum(peaklist$adduct != "", na.rm=TRUE)
cat("  Isotope annotations:", n_isotopes, "\n")
cat("  Adduct annotations:", n_adducts, "\n")

# ============================================================================
# Step 5b: Blank subtraction
# ============================================================================
cat("[Step 5b] Blank subtraction...\n")
blank_idx <- which(meta$sampleType == "blank")
sample_idx <- which(meta$sampleType == "sample")
qc_idx <- which(meta$sampleType == "pool")

# For each feature, compare blank mean to sample mean
blank_means <- rowMeans(fvals[, blank_idx, drop=FALSE], na.rm=TRUE)
sample_means <- rowMeans(fvals[, sample_idx, drop=FALSE], na.rm=TRUE)

# Features where sample signal > 3x blank signal
blank_ratio <- ifelse(blank_means > 0 & !is.na(blank_means), sample_means / blank_means, Inf)
blank_pass <- (blank_ratio > 3 | is.infinite(blank_ratio))
blank_pass[is.na(blank_pass)] <- TRUE  # keep features with no blank data
n_blank_removed <- sum(!blank_pass, na.rm=TRUE)
cat("  Features removed by blank filter:", n_blank_removed, "\n")

# ============================================================================
# Step 5c: QC-based filtering (CV < 30% in QC samples)
# ============================================================================
cat("[Step 5c] QC-based filtering...\n")
qc_vals <- fvals[, qc_idx, drop=FALSE]
qc_cv <- apply(qc_vals, 1, function(x) {
    x <- x[!is.na(x)]
    if(length(x) < 2) return(NA)
    sd(x) / mean(x) * 100
})
qc_pass <- !is.na(qc_cv) & qc_cv < 30
qc_pass[is.na(qc_pass)] <- FALSE
n_qc_removed <- sum(!qc_pass, na.rm=TRUE)
cat("  Features removed by QC CV filter:", n_qc_removed, "\n")

# ============================================================================
# Step 6: Combined filtering — CONVERGE #3
# ============================================================================
cat("[Step 6] Applying combined filters (convergence #3)...\n")
keep <- blank_pass & qc_pass
keep[is.na(keep)] <- FALSE
fvals_filtered <- fvals[keep, ]
fdefs_filtered <- fdefs[keep, ]
n_features_final <- nrow(fvals_filtered)
cat("  Features after all filters:", n_features_final, "\n")

# Save filtered feature table
feature_table <- data.frame(
    feature_id = rownames(fdefs_filtered),
    mz = fdefs_filtered$mzmed,
    rt = fdefs_filtered$rtmed,
    fvals_filtered[, sample_idx, drop=FALSE]
)
write.csv(feature_table, file.path(out_dir, "filtered_features.csv"), row.names=FALSE)

# ============================================================================
# Step 7a: Spectral matching / putative annotation
# ============================================================================
cat("[Step 7a] Putative metabolite annotation by mass matching...\n")
# Use common metabolite masses for negative mode annotation
# Common negative mode adducts: [M-H]-, [M+Cl]-, [M+FA-H]-
# Known urine metabolites with their exact masses
known_metabolites <- data.frame(
    name = c("Hippuric acid", "Citric acid", "Indoxyl sulfate",
             "Creatinine", "Uric acid", "Tryptophan",
             "Phenylacetic acid", "Kynurenic acid", "Xanthurenic acid",
             "Pantothenic acid"),
    mz_neg = c(178.0510, 191.0197, 212.0023,
               112.0510, 167.0211, 203.0826,
               135.0452, 188.0717, 204.0666,
               218.1034),
    stringsAsFactors = FALSE
)

mz_filtered <- fdefs_filtered$mzmed
annotations <- character(n_features_final)
for(i in seq_len(nrow(known_metabolites))) {
    matches <- which(abs(mz_filtered - known_metabolites$mz_neg[i]) < 0.01)
    for(m in matches) {
        annotations[m] <- known_metabolites$name[i]
    }
}
n_annotated <- sum(annotations != "")
cat("  Features with putative annotation:", n_annotated, "\n")

# ============================================================================
# Step 7b: Statistical analysis (fold change, t-test by gender)
# ============================================================================
cat("[Step 7b] Statistical analysis...\n")
# Compare male vs female
male_idx <- which(meta$sampleType == "sample" & meta$gender == "Male")
female_idx <- which(meta$sampleType == "sample" & meta$gender == "Female")

# Map from global indices to sample column positions in fvals_filtered
pvalues <- numeric(n_features_final)
fold_changes <- numeric(n_features_final)
for(i in seq_len(n_features_final)) {
    m_vals <- as.numeric(fvals_filtered[i, male_idx])
    f_vals <- as.numeric(fvals_filtered[i, female_idx])
    m_vals <- m_vals[!is.na(m_vals)]
    f_vals <- f_vals[!is.na(f_vals)]
    if(length(m_vals) >= 2 && length(f_vals) >= 2) {
        tt <- tryCatch(t.test(m_vals, f_vals), error=function(e) NULL)
        if(!is.null(tt)) {
            pvalues[i] <- tt$p.value
            fold_changes[i] <- log2(mean(m_vals) / mean(f_vals))
        } else {
            pvalues[i] <- 1
            fold_changes[i] <- 0
        }
    } else {
        pvalues[i] <- 1
        fold_changes[i] <- 0
    }
}
# FDR correction
padj <- p.adjust(pvalues, method="BH")
n_significant <- sum(padj < 0.05, na.rm=TRUE)
cat("  Significant features (FDR<0.05):", n_significant, "\n")

# ============================================================================
# Step 8: Enrichment/pathway summary — CONVERGE #4
# ============================================================================
cat("[Step 8] Pathway enrichment summary (convergence #4)...\n")
# Combine spectral matching and statistical results
sig_annotated <- sum(padj < 0.05 & annotations != "", na.rm=TRUE)
cat("  Significant AND annotated:", sig_annotated, "\n")

# ============================================================================
# Step 9-10: Generate final report
# ============================================================================
cat("[Step 9-10] Generating report...\n")

# Compute summary metrics
total_scans <- sum(sapply(seq_along(fileNames(raw_data)), function(i) {
    length(scanIndex(filterFile(raw_data, i)))
}))

# mz range
mz_range <- range(fdefs$mzmed)
rt_range <- range(fdefs$rtmed)

# Mean QC CV for final features
mean_qc_cv <- round(mean(qc_cv[keep], na.rm=TRUE), 2)
median_qc_cv <- round(median(qc_cv[keep], na.rm=TRUE), 2)

# Top feature by p-value (regardless of significance)
best_idx <- which.min(pvalues)
top_feature_mz <- round(fdefs_filtered$mzmed[best_idx], 4)
top_feature_pval <- format(pvalues[best_idx], scientific=TRUE, digits=3)
top_feature_fdr <- format(padj[best_idx], scientific=TRUE, digits=3)
top_feature_fc <- round(fold_changes[best_idx], 4)
top_feature_anno <- annotations[best_idx]
if(top_feature_anno == "") top_feature_anno <- "unknown"

# Also use FDR < 0.25 (more permissive for small samples)
n_sig_fdr25 <- sum(padj < 0.25, na.rm=TRUE)

report <- data.frame(
    metric = c(
        "num_samples", "num_qc_pools", "num_blanks",
        "total_scans", "mz_range_min", "mz_range_max",
        "rt_range_min_sec", "rt_range_max_sec",
        "peaks_detected", "features_grouped", "features_after_fill",
        "features_blank_removed", "features_qc_removed",
        "features_final",
        "isotope_annotations", "adduct_annotations",
        "putative_identifications",
        "mean_qc_cv_pct", "median_qc_cv_pct",
        "significant_features_fdr05",
        "significant_features_fdr25",
        "significant_annotated_features",
        "top_feature_mz", "top_feature_pvalue",
        "top_feature_fdr", "top_feature_log2fc",
        "top_feature_annotation"
    ),
    value = c(
        length(sample_idx), length(qc_idx), length(blank_idx),
        total_scans, round(mz_range[1], 4), round(mz_range[2], 4),
        round(rt_range[1], 2), round(rt_range[2], 2),
        nrow(chromPeaks(xdata)), n_features_grouped, n_features_filled,
        n_blank_removed, n_qc_removed,
        n_features_final,
        n_isotopes, n_adducts,
        n_annotated,
        mean_qc_cv, median_qc_cv,
        n_significant,
        n_sig_fdr25,
        sig_annotated,
        top_feature_mz, top_feature_pval,
        top_feature_fdr, top_feature_fc,
        top_feature_anno
    ),
    stringsAsFactors = FALSE
)

write.csv(report, file.path(results_dir, "report.csv"), row.names=FALSE, quote=FALSE)
cat("\nReport written to", file.path(results_dir, "report.csv"), "\n")
cat("\n")
for(i in seq_len(nrow(report))) {
    cat(sprintf("  %s: %s\n", report$metric[i], report$value[i]))
}
REOF

Rscript --no-save "${OUT}/_pipeline.R" "${DATA}" "${OUT}" "${RESULTS}"

fi

echo ""
echo "========================================="
echo "  LC-MS Metabolomics Analysis Complete!"
echo "========================================="
echo ""
cat "${RESULTS}/report.csv"
