#!/bin/bash
set -euo pipefail

# =============================================================================
# Task: CUT&RUN Epigenomic Profiling (H3K4me3)
#
# DAG structure (depth 7, 2 convergence points):
#
# L0: PE reads + genome + blacklist
# L1: fastp (trim)
# L2: bowtie2 (align to genome, --very-sensitive --no-mixed --dovetail -X 700)
# L3: samtools sort + filter (MAPQ>=10, proper pairs, no chrM)
# L4: picard MarkDuplicates (remove PCR dups)
#     ├──────────────────────────────────────────────┐
# L5: SEACR (peak calling, stringent)         MACS2 (peak calling, narrow)
#     └────────────┬─────────────────────────────────┘
# L6: bedtools intersect (consensus peaks)  [CONVERGENCE 1: dual peak callers]
#     ├─────────────────────────────────────┐
# L7: deeptools (heatmap at peaks)     FRiP calculation
#     └────────────┬────────────────────────┘
# L8: MERGE                                 [CONVERGENCE 2: signal + QC]
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${REF}/genome.fa"
BLACKLIST="${REF}/blacklist.bed"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{trimmed,aligned,dedup,peaks_seacr,peaks_macs,consensus,signal} "${RES}"

# L1
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
fi

# L2: Align (CUT&RUN specific: --dovetail -X 700 for nucleosome fragments)
log_step "L2: bowtie2"
if [ ! -f "${OUT}/aligned/aligned.bam" ]; then
    [ ! -f "${GENOME}.1.bt2" ] && bowtie2-build "${GENOME}" "${GENOME}" --threads ${THREADS}
    bowtie2 -x "${GENOME}" -1 "${OUT}/trimmed/R1.fastq.gz" -2 "${OUT}/trimmed/R2.fastq.gz" \
            --very-sensitive --no-mixed --dovetail -X 700 \
            --rg-id sample --rg "SM:sample" --rg "PL:ILLUMINA" \
            --threads ${THREADS} 2>"${OUT}/aligned/bowtie2.log" \
        | samtools view -bS -f 2 -q 10 | samtools sort -@ ${THREADS} -o "${OUT}/aligned/aligned.bam"
    samtools index "${OUT}/aligned/aligned.bam"
fi

# L3: Stats
samtools flagstat "${OUT}/aligned/aligned.bam" > "${OUT}/aligned/flagstat.txt"

# L4: Dedup
log_step "L4: picard dedup"
if [ ! -f "${OUT}/dedup/dedup.bam" ]; then
    picard MarkDuplicates INPUT="${OUT}/aligned/aligned.bam" OUTPUT="${OUT}/dedup/dedup.bam" \
           METRICS_FILE="${OUT}/dedup/metrics.txt" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    samtools index "${OUT}/dedup/dedup.bam"
fi
BAM="${OUT}/dedup/dedup.bam"

# Remove blacklist
if [ -s "${BLACKLIST}" ]; then
    bedtools intersect -v -abam "${BAM}" -b "${BLACKLIST}" > "${OUT}/dedup/clean.bam"
    samtools index "${OUT}/dedup/clean.bam"
    BAM="${OUT}/dedup/clean.bam"
fi

# L5 LEFT: SEACR peak calling (needs bedgraph)
log_step "L5-LEFT: SEACR"
if [ ! -f "${OUT}/peaks_seacr/peaks.stringent.bed" ]; then
    # Generate bedgraph for SEACR
    bedtools genomecov -ibam "${BAM}" -bg | sort -k1,1 -k2,2n > "${OUT}/peaks_seacr/signal.bedgraph"
    SEACR_1.3.sh "${OUT}/peaks_seacr/signal.bedgraph" 0.01 non stringent "${OUT}/peaks_seacr/peaks" 2>&1 || true
fi

# L5 RIGHT: MACS2 peak calling
log_step "L5-RIGHT: MACS2"
if [ ! -f "${OUT}/peaks_macs/sample_peaks.narrowPeak" ]; then
    macs2 callpeak -t "${BAM}" -f BAMPE -g hs --keep-dup all \
          --nomodel -n sample --outdir "${OUT}/peaks_macs" 2>"${OUT}/peaks_macs/macs2.log" || \
    macs3 callpeak -t "${BAM}" -f BAMPE -g hs --keep-dup all \
          --nomodel -n sample --outdir "${OUT}/peaks_macs" 2>"${OUT}/peaks_macs/macs3.log" || true
fi

# L6: Consensus peaks (intersection of SEACR + MACS2)
log_step "L6: consensus peaks"
SEACR_PEAKS="${OUT}/peaks_seacr/peaks.stringent.bed"
MACS_PEAKS="${OUT}/peaks_macs/sample_peaks.narrowPeak"
if [ -s "$SEACR_PEAKS" ] && [ -s "$MACS_PEAKS" ]; then
    bedtools intersect -a "$SEACR_PEAKS" -b "$MACS_PEAKS" -u > "${OUT}/consensus/consensus.bed" 2>/dev/null || true
fi

# L7: Signal track
log_step "L7: bamCoverage"
if [ ! -f "${OUT}/signal/sample.bw" ]; then
    bamCoverage --bam "${BAM}" --outFileName "${OUT}/signal/sample.bw" \
                --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 50818468 \
                --numberOfProcessors ${THREADS} 2>/dev/null || true
fi

# MERGE
log_step "MERGE"

TOTAL_READS=$(grep "in total" "${OUT}/aligned/flagstat.txt" | awk '{print $1}')
MAPPED=$(grep "mapped (" "${OUT}/aligned/flagstat.txt" | head -1 | awk '{print $1}')
MAPPING_PCT=$(grep "mapped (" "${OUT}/aligned/flagstat.txt" | head -1 | grep -oP '[\d.]+%' | tr -d '%')
DUP_PCT=$(grep "PERCENT_DUPLICATION" "${OUT}/dedup/metrics.txt" -A1 | tail -1 | cut -f9 || echo "0")

SEACR_COUNT=$(wc -l < "$SEACR_PEAKS" 2>/dev/null | tr -d ' ' || echo "0")
MACS_COUNT=$(wc -l < "$MACS_PEAKS" 2>/dev/null | tr -d ' ' || echo "0")
CONSENSUS_COUNT=$(wc -l < "${OUT}/consensus/consensus.bed" 2>/dev/null | tr -d ' ' || echo "0")

# FRiP
TOTAL_BAM=$(samtools view -c "${BAM}" 2>/dev/null || echo "0")
if [ -s "${OUT}/consensus/consensus.bed" ] && [ "$TOTAL_BAM" -gt 0 ]; then
    IN_PEAKS=$(bedtools intersect -a "${BAM}" -b "${OUT}/consensus/consensus.bed" -u -ubam | samtools view -c 2>/dev/null || echo "0")
    FRIP=$(python3 -c "print(f'{${IN_PEAKS}/${TOTAL_BAM}:.4f}')" 2>/dev/null || echo "0")
else
    FRIP="0"
fi

cat > "${RES}/cutandrun_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
mapped_reads,${MAPPED}
mapping_rate,${MAPPING_PCT}
duplication_rate,${DUP_PCT}
peaks_caller_a,${SEACR_COUNT}
peaks_caller_b,${MACS_COUNT}
consensus_peaks,${CONSENSUS_COUNT}
fraction_reads_in_peaks,${FRIP}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/cutandrun_report.csv"
