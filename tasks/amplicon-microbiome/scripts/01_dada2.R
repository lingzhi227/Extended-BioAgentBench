library(dada2)

samples <- c("sample1", "sample1a", "sample2", "sample2a")
trimdir <- "outputs/trimmed"

fnFs <- file.path(trimdir, paste0(samples, "_R1.fastq.gz"))
fnRs <- file.path(trimdir, paste0(samples, "_R2.fastq.gz"))
names(fnFs) <- samples; names(fnRs) <- samples

filtFs <- file.path("outputs/filtered", paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path("outputs/filtered", paste0(samples, "_R_filt.fastq.gz"))
names(filtFs) <- samples; names(filtRs) <- samples

# Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(230, 200), maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
cat("Filter:\n"); print(out)
write.csv(out, "outputs/dada2/filter_stats.csv")

keep <- out[,"reads.out"] > 0
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]

# Learn errors
cat("Learning errors...\n")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Denoise
cat("Denoising...\n")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge
cat("Merging...\n")
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Sequence table + chimera removal
seqtab <- makeSequenceTable(merged)
cat("Seq table:", dim(seqtab), "\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
cat("Non-chimeric ASVs:", ncol(seqtab.nochim), "\n")

# Track reads
track <- cbind(out[keep,],
               sapply(dadaFs, function(x) sum(getUniques(x))),
               sapply(merged, function(x) sum(getUniques(x))),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nonchim")
write.csv(track, "outputs/dada2/track_reads.csv")

# Write ASV sequences as FASTA
seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("ASV", seq_along(seqs))
fasta <- character(2 * length(seqs))
for (i in seq_along(seqs)) {
  fasta[2*i - 1] <- paste0(">", asv_ids[i])
  fasta[2*i] <- seqs[i]
}
writeLines(fasta, "outputs/dada2/asv_seqs.fasta")

# Write ASV count table (samples x ASVs)
colnames(seqtab.nochim) <- asv_ids
df <- as.data.frame(seqtab.nochim)
df$sample_id <- rownames(df)
df <- df[, c("sample_id", asv_ids)]
write.table(df, "outputs/dada2/asv_table.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# Summary
cat("total_input_reads:", sum(out[,"reads.in"]), "\n", file="outputs/dada2/summary.txt")
cat("total_filtered_reads:", sum(out[keep,"reads.out"]), "\n", file="outputs/dada2/summary.txt", append=TRUE)
cat("total_asvs:", ncol(seqtab.nochim), "\n", file="outputs/dada2/summary.txt", append=TRUE)
cat("total_chimeras:", ncol(seqtab) - ncol(seqtab.nochim), "\n", file="outputs/dada2/summary.txt", append=TRUE)
cat("Done!\n")
