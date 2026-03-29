library(dada2)

# Read ASV sequences
lines <- readLines("outputs/dada2/asv_seqs.fasta")
ids <- lines[seq(1, length(lines), 2)]
seqs <- lines[seq(2, length(lines), 2)]
names(seqs) <- sub("^>", "", ids)

cat("Assigning taxonomy to", length(seqs), "ASVs...\n")

# Assign taxonomy to genus level
taxa <- assignTaxonomy(seqs,
                       "reference/silva_nr99_v138.1_train_set.fa.gz",
                       multithread = TRUE, verbose = TRUE)

# Create output (genus-level classification)
tax_df <- as.data.frame(taxa)
tax_df$ASV_ID <- names(seqs)
tax_df <- tax_df[, c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
write.table(tax_df, "outputs/taxonomy/taxonomy.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# Summary
cat("Classified to phylum:", sum(!is.na(tax_df$Phylum)), "of", nrow(tax_df), "\n")
cat("Classified to genus:", sum(!is.na(tax_df$Genus)), "of", nrow(tax_df), "\n")

# Phylum summary
phylum_tab <- table(tax_df$Phylum, useNA="no")
phylum_df <- data.frame(Phylum=names(phylum_tab), Count=as.integer(phylum_tab))
phylum_df <- phylum_df[order(-phylum_df$Count),]
write.table(phylum_df, "outputs/taxonomy/phylum_summary.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("Top phyla:\n")
print(head(phylum_df, 10))
