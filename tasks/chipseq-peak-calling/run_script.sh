#!/bin/bash
set -e

# =============================================================================
# ChIP-seq Peak Calling: TAL1 Binding Site Identification
# Pipeline: trimmomatic -> bwa -> samtools -> macs2 -> bedtools -> deeptools
# Data: Mouse G1E cells + Megakaryocytes, TAL1 ChIP-seq (chr19 subset)
# Source: Galaxy Training Network / Wu et al. 2014 (GEO GSE51338)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
REF_DIR="${SCRIPT_DIR}/reference"
OUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "=================================================================="
    echo "Started at: $(date)"
}

mkdir -p "${OUT_DIR}/trimmed" "${OUT_DIR}/aligned" "${OUT_DIR}/peaks" "${OUT_DIR}/comparison"
mkdir -p "${RESULTS_DIR}"

REF="${REF_DIR}/mm10_chr19.fa"

# ==========================================================================
# STEP 1: Index reference genome (bwa)
# ==========================================================================
log_step "Indexing reference with bwa"
bwa index "${REF}"
samtools faidx "${REF}"

# ==========================================================================
# STEP 2: Trim reads (trimmomatic)
# ==========================================================================
log_step "Trimming reads with trimmomatic"

for SAMPLE in G1E_input_R1 G1E_input_R2 G1E_Tal1_R1 G1E_Tal1_R2 \
              Mega_input_R1 Mega_input_R2 Mega_Tal1_R1 Mega_Tal1_R2; do
    echo "  Trimming ${SAMPLE}..."
    trimmomatic SE -threads ${THREADS} \
        "${DATA_DIR}/${SAMPLE}.fastq" \
        "${OUT_DIR}/trimmed/${SAMPLE}.trimmed.fastq" \
        SLIDINGWINDOW:4:20
done

# ==========================================================================
# STEP 3: Align reads to reference (bwa mem)
# ==========================================================================
log_step "Aligning reads with bwa mem"

for SAMPLE in G1E_input_R1 G1E_input_R2 G1E_Tal1_R1 G1E_Tal1_R2 \
              Mega_input_R1 Mega_input_R2 Mega_Tal1_R1 Mega_Tal1_R2; do
    echo "  Aligning ${SAMPLE}..."
    bwa mem -t ${THREADS} "${REF}" "${OUT_DIR}/trimmed/${SAMPLE}.trimmed.fastq" \
        | samtools sort -@ ${THREADS} -o "${OUT_DIR}/aligned/${SAMPLE}.bam"
    samtools index "${OUT_DIR}/aligned/${SAMPLE}.bam"
done

# ==========================================================================
# STEP 4: Alignment statistics (samtools idxstats)
# ==========================================================================
log_step "Generating alignment statistics with samtools"

echo -e "sample\tref\tlength\tmapped\tunmapped" > "${RESULTS_DIR}/alignment_stats.tsv"
for SAMPLE in G1E_input_R1 G1E_input_R2 G1E_Tal1_R1 G1E_Tal1_R2 \
              Mega_input_R1 Mega_input_R2 Mega_Tal1_R1 Mega_Tal1_R2; do
    samtools idxstats "${OUT_DIR}/aligned/${SAMPLE}.bam" \
        | awk -v s="${SAMPLE}" 'BEGIN{OFS="\t"} {print s, $1, $2, $3, $4}' \
        >> "${RESULTS_DIR}/alignment_stats.tsv"
done

# ==========================================================================
# STEP 5: Peak calling with MACS2 — G1E cells
# ==========================================================================
log_step "Calling peaks for G1E with macs2"

macs2 callpeak \
    -t "${OUT_DIR}/aligned/G1E_Tal1_R1.bam" "${OUT_DIR}/aligned/G1E_Tal1_R2.bam" \
    -c "${OUT_DIR}/aligned/G1E_input_R1.bam" "${OUT_DIR}/aligned/G1E_input_R2.bam" \
    -f BAM -g mm --call-summits \
    -n G1E_TAL1 --outdir "${OUT_DIR}/peaks/"

# ==========================================================================
# STEP 6: Peak calling with MACS2 — Megakaryocytes
# ==========================================================================
log_step "Calling peaks for Megakaryocytes with macs2"

macs2 callpeak \
    -t "${OUT_DIR}/aligned/Mega_Tal1_R1.bam" "${OUT_DIR}/aligned/Mega_Tal1_R2.bam" \
    -c "${OUT_DIR}/aligned/Mega_input_R1.bam" "${OUT_DIR}/aligned/Mega_input_R2.bam" \
    -f BAM -g mm --call-summits \
    -n Mega_TAL1 --outdir "${OUT_DIR}/peaks/"

# ==========================================================================
# STEP 7: Compare peaks between cell types (bedtools intersect)
# ==========================================================================
log_step "Comparing peaks with bedtools"

# Common peaks (shared between G1E and Megakaryocytes)
bedtools intersect \
    -a "${OUT_DIR}/peaks/G1E_TAL1_peaks.narrowPeak" \
    -b "${OUT_DIR}/peaks/Mega_TAL1_peaks.narrowPeak" \
    > "${OUT_DIR}/comparison/common_peaks.bed"

# G1E-unique peaks
bedtools intersect \
    -a "${OUT_DIR}/peaks/G1E_TAL1_peaks.narrowPeak" \
    -b "${OUT_DIR}/peaks/Mega_TAL1_peaks.narrowPeak" \
    -v > "${OUT_DIR}/comparison/g1e_unique_peaks.bed"

# Megakaryocyte-unique peaks
bedtools intersect \
    -a "${OUT_DIR}/peaks/Mega_TAL1_peaks.narrowPeak" \
    -b "${OUT_DIR}/peaks/G1E_TAL1_peaks.narrowPeak" \
    -v > "${OUT_DIR}/comparison/mega_unique_peaks.bed"

# ==========================================================================
# STEP 8: Normalized signal tracks (deeptools bamCompare)
# ==========================================================================
log_step "Computing normalized signal with deeptools bamCompare"

bamCompare \
    -b1 "${OUT_DIR}/aligned/G1E_Tal1_R1.bam" \
    -b2 "${OUT_DIR}/aligned/G1E_input_R1.bam" \
    --operation log2 --binSize 50 \
    -p ${THREADS} \
    -o "${OUT_DIR}/comparison/G1E_R1_log2ratio.bw"

bamCompare \
    -b1 "${OUT_DIR}/aligned/Mega_Tal1_R1.bam" \
    -b2 "${OUT_DIR}/aligned/Mega_input_R1.bam" \
    --operation log2 --binSize 50 \
    -p ${THREADS} \
    -o "${OUT_DIR}/comparison/Mega_R1_log2ratio.bw"

# ==========================================================================
# STEP 9: Generate final results CSV
# ==========================================================================
log_step "Generating final results"

# Main output: peak comparison summary with coordinates and scores
echo "chrom,start,end,name,score,strand,signal_value,pvalue,qvalue,peak,cell_type,status" \
    > "${RESULTS_DIR}/peak_comparison.csv"

# G1E peaks — mark shared vs unique
awk 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"G1E","shared"}' \
    "${OUT_DIR}/comparison/common_peaks.bed" >> "${RESULTS_DIR}/peak_comparison.csv"
awk 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"G1E","unique"}' \
    "${OUT_DIR}/comparison/g1e_unique_peaks.bed" >> "${RESULTS_DIR}/peak_comparison.csv"

# Megakaryocyte peaks — shared vs unique
awk 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"Megakaryocyte","shared"}' \
    "${OUT_DIR}/comparison/common_peaks.bed" >> "${RESULTS_DIR}/peak_comparison.csv"
awk 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"Megakaryocyte","unique"}' \
    "${OUT_DIR}/comparison/mega_unique_peaks.bed" >> "${RESULTS_DIR}/peak_comparison.csv"

# Summary stats
TOTAL_G1E=$(wc -l < "${OUT_DIR}/peaks/G1E_TAL1_peaks.narrowPeak")
TOTAL_MEGA=$(wc -l < "${OUT_DIR}/peaks/Mega_TAL1_peaks.narrowPeak")
COMMON=$(wc -l < "${OUT_DIR}/comparison/common_peaks.bed")
G1E_UNIQUE=$(wc -l < "${OUT_DIR}/comparison/g1e_unique_peaks.bed")
MEGA_UNIQUE=$(wc -l < "${OUT_DIR}/comparison/mega_unique_peaks.bed")

echo "metric,value" > "${RESULTS_DIR}/peak_summary.csv"
echo "g1e_total_peaks,${TOTAL_G1E}" >> "${RESULTS_DIR}/peak_summary.csv"
echo "mega_total_peaks,${TOTAL_MEGA}" >> "${RESULTS_DIR}/peak_summary.csv"
echo "common_peaks,${COMMON}" >> "${RESULTS_DIR}/peak_summary.csv"
echo "g1e_unique_peaks,${G1E_UNIQUE}" >> "${RESULTS_DIR}/peak_summary.csv"
echo "mega_unique_peaks,${MEGA_UNIQUE}" >> "${RESULTS_DIR}/peak_summary.csv"

echo ""
echo "=== Pipeline complete ==="
echo "Results in: ${RESULTS_DIR}/"
ls -lh "${RESULTS_DIR}/"
cat "${RESULTS_DIR}/peak_summary.csv"
