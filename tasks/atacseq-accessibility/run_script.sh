#!/bin/bash
set -e

# =============================================================================
# Task 15: ATAC-seq Chromatin Accessibility Analysis (depth 8)
#
# L0: reads
# L1: trim_galore (adapter trim)
# L2: bowtie2 (align)
# L3: samtools (sort, index, filter chrM, flagstat)
# L4: picard MarkDuplicates (dedup)
# L5: bedtools intersect -v (blacklist removal)
# L6: deeptools alignmentSieve --ATACshift (Tn5 correction)
#     ├──────────────────────────────┐
# L7: macs2 callpeak (BAMPE)      deeptools bamCoverage (bigWig)
#     ├──────────┐
# L8: homer      featureCounts (FRiP)
#     findMotifs
#     MERGE (QC report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${REF}/chr22.fa"
BLACKLIST="${REF}/blacklist_chr22.bed"
GENES="${REF}/genes_chr22.gtf"
SAMPLE="SRR891268"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,aligned,dedup,filtered,shifted,peaks,signal,motifs,frip,qc} "${RES}"

# --- Build bowtie2 index ---
log_step "INDEX: bowtie2-build"
if [ ! -f "${REF}/chr22.1.bt2" ]; then
    bowtie2-build "${GENOME}" "${REF}/chr22" --threads ${THREADS}
    samtools faidx "${GENOME}"
    cut -f1,2 "${GENOME}.fai" > "${REF}/chr22.chrom.sizes"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L1: Adapter trimming with Trim Galore
# ===========================================================================
log_step "L1: Trim Galore adapter trimming"
if [ ! -f "${OUT}/trimmed/${SAMPLE}_R1_val_1.fq.gz" ]; then
    trim_galore --paired --fastqc --cores ${THREADS} \
                -o "${OUT}/trimmed" \
                "${DATA}/${SAMPLE}_R1.fastq.gz" "${DATA}/${SAMPLE}_R2.fastq.gz"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2: Alignment with Bowtie2
# ===========================================================================
log_step "L2: Bowtie2 alignment"
if [ ! -f "${OUT}/aligned/${SAMPLE}.bam" ]; then
    bowtie2 -x "${REF}/chr22" \
            -1 "${OUT}/trimmed/${SAMPLE}_R1_val_1.fq.gz" \
            -2 "${OUT}/trimmed/${SAMPLE}_R2_val_2.fq.gz" \
            --very-sensitive -X 2000 --no-mixed --no-discordant \
            --rg-id "${SAMPLE}" --rg "SM:${SAMPLE}" --rg "PL:ILLUMINA" \
            --threads ${THREADS} 2>"${OUT}/aligned/bowtie2.log" \
        | samtools view -bS -f 2 -q 30 - > "${OUT}/aligned/${SAMPLE}.bam"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3: Sort, index, filter mitochondrial reads, flagstat
# ===========================================================================
log_step "L3: samtools sort/index/filter"
if [ ! -f "${OUT}/aligned/${SAMPLE}.sorted.bam" ]; then
    samtools sort -@ ${THREADS} -o "${OUT}/aligned/${SAMPLE}.sorted.bam" "${OUT}/aligned/${SAMPLE}.bam"
    samtools index "${OUT}/aligned/${SAMPLE}.sorted.bam"
    samtools flagstat "${OUT}/aligned/${SAMPLE}.sorted.bam" > "${OUT}/qc/flagstat.txt"
    # Filter chrM (not applicable for chr22-only data, but shown for correctness)
    samtools idxstats "${OUT}/aligned/${SAMPLE}.sorted.bam" > "${OUT}/qc/idxstats.txt"
    CHRM_READS=$(awk '$1=="chrM" {print $3}' "${OUT}/qc/idxstats.txt" || echo "0")
    echo "chrM reads: ${CHRM_READS}" | tee "${OUT}/qc/chrM_filter.txt"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4: Mark and remove PCR duplicates with Picard
# ===========================================================================
log_step "L4: Picard MarkDuplicates"
if [ ! -f "${OUT}/dedup/${SAMPLE}.dedup.bam" ]; then
    picard MarkDuplicates \
        INPUT="${OUT}/aligned/${SAMPLE}.sorted.bam" \
        OUTPUT="${OUT}/dedup/${SAMPLE}.dedup.bam" \
        METRICS_FILE="${OUT}/qc/picard_dedup_metrics.txt" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT
    samtools index "${OUT}/dedup/${SAMPLE}.dedup.bam"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5: Remove ENCODE blacklist regions
# ===========================================================================
log_step "L5: Blacklist region removal"
if [ ! -f "${OUT}/filtered/${SAMPLE}.filtered.bam" ]; then
    if [ -s "${BLACKLIST}" ]; then
        bedtools intersect -v -abam "${OUT}/dedup/${SAMPLE}.dedup.bam" \
                           -b "${BLACKLIST}" > "${OUT}/filtered/${SAMPLE}.filtered.bam"
    else
        cp "${OUT}/dedup/${SAMPLE}.dedup.bam" "${OUT}/filtered/${SAMPLE}.filtered.bam"
    fi
    samtools index "${OUT}/filtered/${SAMPLE}.filtered.bam"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L6: Tn5 transposase shift correction (+4/-5 bp)
# ===========================================================================
log_step "L6: Tn5 shift correction (deeptools alignmentSieve)"
if [ ! -f "${OUT}/shifted/${SAMPLE}.shifted.bam" ]; then
    alignmentSieve --bam "${OUT}/filtered/${SAMPLE}.filtered.bam" \
                   --outFile "${OUT}/shifted/${SAMPLE}.shifted.bam" \
                   --ATACshift \
                   --numberOfProcessors ${THREADS}
    samtools sort -@ ${THREADS} -o "${OUT}/shifted/${SAMPLE}.shifted.sorted.bam" \
                  "${OUT}/shifted/${SAMPLE}.shifted.bam"
    samtools index "${OUT}/shifted/${SAMPLE}.shifted.sorted.bam"
else echo "Skipping (exists)"; fi
FINAL_BAM="${OUT}/shifted/${SAMPLE}.shifted.sorted.bam"

# --- Fragment size distribution (ATAC-specific QC) ---
log_step "L6-QC: Fragment size distribution"
if [ ! -f "${OUT}/qc/fragment_sizes.png" ]; then
    bamPEFragmentSizes --bamfiles "${FINAL_BAM}" \
                       --histogram "${OUT}/qc/fragment_sizes.png" \
                       --outRawFragmentLengths "${OUT}/qc/fragment_lengths.txt" \
                       --numberOfProcessors ${THREADS} 2>/dev/null || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L7 LEFT: Peak calling with MACS2 (BAMPE mode)
# ===========================================================================
log_step "L7-LEFT: MACS2 peak calling (BAMPE)"
if [ ! -f "${OUT}/peaks/${SAMPLE}_peaks.narrowPeak" ]; then
    macs3 callpeak -t "${FINAL_BAM}" \
                   -f BAMPE -g hs --keep-dup all \
                   --nomodel --shift -100 --extsize 200 \
                   --call-summits \
                   -n "${SAMPLE}" --outdir "${OUT}/peaks" \
                   2>"${OUT}/peaks/macs2.log"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L7 RIGHT: Signal track with deeptools bamCoverage
# ===========================================================================
log_step "L7-RIGHT: deeptools bamCoverage (bigWig)"
if [ ! -f "${OUT}/signal/${SAMPLE}.bw" ]; then
    bamCoverage --bam "${FINAL_BAM}" \
                --outFileName "${OUT}/signal/${SAMPLE}.bw" \
                --binSize 10 --normalizeUsing RPGC \
                --effectiveGenomeSize 50818468 \
                --numberOfProcessors ${THREADS}
else echo "Skipping (exists)"; fi

# ===========================================================================
# L8 LEFT: Motif enrichment with HOMER
# ===========================================================================
log_step "L8-LEFT: HOMER motif enrichment"
if [ ! -f "${OUT}/motifs/knownResults.html" ]; then
    findMotifsGenome.pl "${OUT}/peaks/${SAMPLE}_summits.bed" \
                        "${GENOME}" "${OUT}/motifs" \
                        -size 200 -mask -p ${THREADS} 2>&1 || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L8 RIGHT: FRiP (Fraction of Reads in Peaks)
# ===========================================================================
log_step "L8-RIGHT: FRiP calculation"
TOTAL_READS=$(samtools view -c "${FINAL_BAM}" 2>/dev/null || echo "0")
if [ -f "${OUT}/peaks/${SAMPLE}_peaks.narrowPeak" ]; then
    # Convert peaks to SAF for featureCounts
    awk 'BEGIN{OFS="\t"} {print $4,$1,$2,$3,"."}' \
        "${OUT}/peaks/${SAMPLE}_peaks.narrowPeak" > "${OUT}/frip/peaks.saf"
    featureCounts -a "${OUT}/frip/peaks.saf" -F SAF \
                  -p --countReadPairs \
                  -o "${OUT}/frip/featureCounts.txt" \
                  "${FINAL_BAM}" 2>&1 || true
    READS_IN_PEAKS=$(tail -n +3 "${OUT}/frip/featureCounts.txt" 2>/dev/null \
        | awk '{sum+=$NF} END{print sum}' || echo "0")
    if [ "${TOTAL_READS}" -gt 0 ] 2>/dev/null; then
        FRIP=$(python3 -c "print(f'{${READS_IN_PEAKS}/${TOTAL_READS}:.4f}')")
    else
        FRIP="N/A"
    fi
else
    READS_IN_PEAKS=0
    FRIP="N/A"
fi
echo "FRiP: ${FRIP} (${READS_IN_PEAKS}/${TOTAL_READS})"

# ===========================================================================
# L8: TSS enrichment score
# ===========================================================================
log_step "L8: TSS enrichment"
if [ ! -f "${OUT}/qc/tss_enrichment.png" ]; then
    computeMatrix reference-point -S "${OUT}/signal/${SAMPLE}.bw" \
                                  -R "${GENES}" \
                                  --referencePoint TSS \
                                  -b 2000 -a 2000 \
                                  -o "${OUT}/qc/tss_matrix.gz" \
                                  --numberOfProcessors ${THREADS} 2>/dev/null || true
    plotProfile -m "${OUT}/qc/tss_matrix.gz" \
                -out "${OUT}/qc/tss_enrichment.png" \
                --outFileNameData "${OUT}/qc/tss_enrichment.tab" 2>/dev/null || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# MERGE: Comprehensive QC report
# ===========================================================================
log_step "MERGE: Building results"

NUM_PEAKS=$(wc -l < "${OUT}/peaks/${SAMPLE}_peaks.narrowPeak" 2>/dev/null | tr -d ' ' || echo "0")
DEDUP_PCT=$(grep "PERCENT_DUPLICATION" "${OUT}/qc/picard_dedup_metrics.txt" -A1 | tail -1 | cut -f9 || echo "N/A")
MAPPED_READS=$(grep "mapped (" "${OUT}/qc/flagstat.txt" | head -1 | awk '{print $1}' || echo "N/A")
ALIGN_RATE=$(grep "overall alignment rate" "${OUT}/aligned/bowtie2.log" | awk '{print $1}' || echo "N/A")
NUM_MOTIFS=$(ls "${OUT}/motifs/knownResults/"*.motif 2>/dev/null | wc -l | tr -d ' ' || echo "0")

cat > "${RES}/atacseq_report.csv" << CSVEOF
metric,value
sample,${SAMPLE}
total_reads,$(( TOTAL_READS ))
alignment_rate,${ALIGN_RATE}
mapped_reads,${MAPPED_READS}
duplication_rate,${DEDUP_PCT}
chrM_reads,${CHRM_READS}
num_peaks,${NUM_PEAKS}
frip,${FRIP}
known_motifs_enriched,${NUM_MOTIFS}
CSVEOF

# Copy peak file to results
cp "${OUT}/peaks/${SAMPLE}_peaks.narrowPeak" "${RES}/peaks.narrowPeak" 2>/dev/null || true

# Top 10 motifs
if [ -f "${OUT}/motifs/knownResults.txt" ]; then
    head -11 "${OUT}/motifs/knownResults.txt" > "${RES}/top_motifs.tsv"
fi

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/atacseq_report.csv"
echo ""
ls -lh "${RES}/"
