# Compile final report from all pipeline outputs

# Read DADA2 summary
dada2_summary <- readLines("outputs/dada2/summary.txt")
parse_val <- function(lines, key) {
  line <- grep(paste0("^", key, ":"), lines, value=TRUE)
  trimws(sub(paste0("^", key, ": *"), "", line))
}

total_input <- parse_val(dada2_summary, "total_input_reads")
total_filtered <- parse_val(dada2_summary, "total_filtered_reads")
total_asvs <- parse_val(dada2_summary, "total_asvs")
total_chimeras <- parse_val(dada2_summary, "total_chimeras")

# Taxonomy
tax <- read.table("outputs/taxonomy/taxonomy.tsv", header=TRUE, sep="\t", fill=TRUE)
classified <- sum(!is.na(tax$Phylum))
phylum_tab <- table(tax$Phylum, useNA="no")
top_phylum <- names(sort(phylum_tab, decreasing=TRUE))[1]
top_phylum_pct <- round(100 * max(phylum_tab) / sum(phylum_tab), 1)

# Diversity
div <- read.table("outputs/diversity/diversity_metrics.tsv", header=TRUE, sep="\t")
mean_shannon <- div$value[div$metric == "mean_shannon"]
mean_simpson <- div$value[div$metric == "mean_simpson"]
mean_observed <- div$value[div$metric == "mean_observed_richness"]
mean_pd <- div$value[div$metric == "mean_faith_pd"]

# DA results
da <- read.table("outputs/diversity/da_results.tsv", header=TRUE, sep="\t")
da_features <- sum(da$padj < 0.05, na.rm=TRUE)

# PICRUSt2 pathways
pw_count <- 0
pw_files <- c(
  "outputs/picrust2/out/pathways_out/path_abun_unstrat.tsv.gz",
  "outputs/picrust2/out/pathways_out/path_abun_unstrat.tsv"
)
for (pf in pw_files) {
  if (file.exists(pf)) {
    if (grepl("\\.gz$", pf)) {
      pw <- read.table(gzfile(pf), header=TRUE, sep="\t", check.names=FALSE)
    } else {
      pw <- read.table(pf, header=TRUE, sep="\t", check.names=FALSE)
    }
    pw_count <- nrow(pw)
    break
  }
}

# Write report
report <- data.frame(
  metric = c("total_samples", "total_input_reads", "total_filtered_reads",
             "total_sequence_variants", "chimeras_removed",
             "classified_variants", "top_phylum", "top_phylum_pct",
             "mean_shannon", "mean_simpson", "mean_observed_richness",
             "mean_phylogenetic_diversity", "predicted_pathways",
             "differentially_abundant_features"),
  value = c(4, total_input, total_filtered, total_asvs, total_chimeras,
            classified, top_phylum, top_phylum_pct,
            mean_shannon, mean_simpson, mean_observed, mean_pd,
            pw_count, da_features)
)
write.csv(report, "results/report.csv", row.names=FALSE, quote=FALSE)
cat("=== Final Report ===\n")
for (i in seq_len(nrow(report))) {
  cat(sprintf("  %s = %s\n", report$metric[i], report$value[i]))
}
