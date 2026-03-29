library(phyloseq)

# Load data
asv_tab <- read.table("outputs/dada2/asv_table.tsv", header=TRUE, sep="\t",
                       row.names=1, check.names=FALSE)
asv_mat <- as.matrix(asv_tab)
tax_tab <- read.table("outputs/taxonomy/taxonomy.tsv", header=TRUE, sep="\t",
                       row.names=1, check.names=FALSE, fill=TRUE)
tax_mat <- as.matrix(tax_tab)
meta <- read.table("data/metadata.tsv", header=TRUE, sep="\t", row.names=1)

ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows=FALSE),
  tax_table(tax_mat),
  sample_data(meta)
)

# Differential abundance using Wilcoxon test per ASV (relative abundance)
otu <- as(otu_table(ps), "matrix")
groups <- sample_data(ps)$treatment
taxa <- tax_table(ps)

# Convert to relative abundance
rel_abund <- sweep(otu, 1, rowSums(otu), "/")

results <- data.frame(ASV=character(), pvalue=numeric(), log2fc=numeric(),
                      Phylum=character(), stringsAsFactors=FALSE)

for (i in seq_len(ncol(rel_abund))) {
  grp1 <- rel_abund[groups == "control", i]
  grp2 <- rel_abund[groups == "treated", i]
  if (sum(c(grp1, grp2) > 0) >= 2) {
    if (var(c(grp1, grp2)) > 0) {
      wt <- suppressWarnings(wilcox.test(grp1, grp2))
      fc <- log2((mean(grp2) + 1e-6) / (mean(grp1) + 1e-6))
      results <- rbind(results, data.frame(
        ASV = colnames(rel_abund)[i],
        pvalue = wt$p.value,
        log2fc = round(fc, 4),
        Phylum = as.character(taxa[colnames(rel_abund)[i], "Phylum"])
      ))
    }
  }
}

# Adjust p-values
if (nrow(results) > 0) {
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$pvalue),]
}

write.table(results, "outputs/diversity/da_results.tsv", sep="\t",
            row.names=FALSE, quote=FALSE)
cat("Tested ASVs:", nrow(results), "\n")
cat("DA features (padj<0.05):", sum(results$padj < 0.05, na.rm=TRUE), "\n")
cat("DA features (p<0.05):", sum(results$pvalue < 0.05, na.rm=TRUE), "\n")
