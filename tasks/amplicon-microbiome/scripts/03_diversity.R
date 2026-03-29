library(phyloseq)
library(vegan)
library(ape)

# Load ASV table
asv_tab <- read.table("outputs/dada2/asv_table.tsv", header=TRUE, sep="\t",
                       row.names=1, check.names=FALSE)
asv_mat <- as.matrix(asv_tab)

# Load taxonomy
tax_tab <- read.table("outputs/taxonomy/taxonomy.tsv", header=TRUE, sep="\t",
                       row.names=1, check.names=FALSE, fill=TRUE)
tax_mat <- as.matrix(tax_tab)

# Load tree
tree <- read.tree("outputs/phylogeny/tree.nwk")

# Load metadata
meta <- read.table("data/metadata.tsv", header=TRUE, sep="\t", row.names=1)

# Build phyloseq object
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows=FALSE),
  tax_table(tax_mat),
  sample_data(meta),
  phy_tree(tree)
)
cat("Phyloseq object:\n")
print(ps)

# Alpha diversity
alpha <- estimate_richness(ps, measures=c("Observed", "Shannon", "Simpson"))
alpha$sample_id <- rownames(alpha)

# Faith's phylogenetic diversity
faith_pd <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  tree <- phy_tree(ps)
  pd_vals <- numeric(nrow(otu))
  for (i in seq_len(nrow(otu))) {
    present <- colnames(otu)[otu[i,] > 0]
    if (length(present) > 1) {
      subtree <- keep.tip(tree, present)
      pd_vals[i] <- sum(subtree$edge.length)
    } else {
      pd_vals[i] <- 0
    }
  }
  names(pd_vals) <- rownames(otu)
  return(pd_vals)
}

pd <- faith_pd(ps)
alpha$faith_pd <- pd[rownames(alpha)]

cat("Alpha diversity:\n")
print(alpha)
write.table(alpha, "outputs/diversity/alpha_diversity.tsv", sep="\t",
            row.names=FALSE, quote=FALSE)

# Beta diversity: weighted UniFrac
cat("Computing UniFrac...\n")
unifrac_dist <- UniFrac(ps, weighted=TRUE)
write.table(as.matrix(unifrac_dist), "outputs/diversity/weighted_unifrac.tsv",
            sep="\t", quote=FALSE)

# Bray-Curtis
bray_dist <- vegdist(asv_mat, method="bray")
write.table(as.matrix(bray_dist), "outputs/diversity/bray_curtis.tsv",
            sep="\t", quote=FALSE)

# Summary
div_summary <- data.frame(
  metric = c("mean_shannon", "mean_simpson", "mean_observed_richness",
             "mean_faith_pd"),
  value = c(round(mean(alpha$Shannon), 4),
            round(mean(alpha$Simpson), 4),
            round(mean(alpha$Observed), 0),
            round(mean(alpha$faith_pd), 4))
)
write.table(div_summary, "outputs/diversity/diversity_metrics.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
cat("Diversity summary:\n")
print(div_summary)
