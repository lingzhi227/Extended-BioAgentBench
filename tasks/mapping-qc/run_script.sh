#!/bin/bash
set -e

# =============================================================================
# Task 23: Genome Coverage and Mapping Quality Analysis
#
# DAG (depth 5):
#
# L0: PE reads + reference
# L1: fastp (trim)
# L2: bwa mem (align)
#     ├──────────────────────────────────────┐
# L3: samtools stats/flagstat      mosdepth (coverage)
#     │                                │
# L4: samtools idxstats          ├── per-base coverage
#     (per-contig mapping)       └── per-window coverage
#     │                                │
#     └── bedtools genomecov     deeptools bamCoverage
#         (coverage histogram)    (normalized signal)
# L5: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

REFERENCE="${REF}/reference.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,aligned,stats,coverage} "${RES}"

# ===========================================================================
# L1: Trimming
# ===========================================================================
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq" --in2 "${DATA}/reads_R2.fastq" \
          --out1 "${OUT}/trimmed/R1.fastq" --out2 "${OUT}/trimmed/R2.fastq" \
          --detect_adapter_for_pe --thread ${THREADS} \
          --json "${OUT}/trimmed/fastp.json"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2: Alignment
# ===========================================================================
log_step "L2: bwa mem"
if [ ! -f "${OUT}/aligned/reads.sorted.bam" ]; then
    bwa index "${REFERENCE}" 2>/dev/null
    bwa mem -t ${THREADS} -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
            "${REFERENCE}" "${OUT}/trimmed/R1.fastq" "${OUT}/trimmed/R2.fastq" \
        | samtools sort -@ ${THREADS} -o "${OUT}/aligned/reads.sorted.bam"
    samtools index "${OUT}/aligned/reads.sorted.bam"
else echo "Skipping (exists)"; fi
BAM="${OUT}/aligned/reads.sorted.bam"

# ===========================================================================
# L3-A: samtools stats + flagstat
# ===========================================================================
log_step "L3-A: samtools stats"
if [ ! -f "${OUT}/stats/samtools_stats.txt" ]; then
    samtools stats "${BAM}" > "${OUT}/stats/samtools_stats.txt"
    samtools flagstat "${BAM}" > "${OUT}/stats/flagstat.txt"
    samtools idxstats "${BAM}" > "${OUT}/stats/idxstats.txt"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3-B: mosdepth coverage analysis
# ===========================================================================
log_step "L3-B: mosdepth"
if [ ! -f "${OUT}/coverage/sample.mosdepth.summary.txt" ]; then
    mosdepth -t ${THREADS} --by 1000 "${OUT}/coverage/sample" "${BAM}"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4-A: bedtools genomecov (coverage histogram)
# ===========================================================================
log_step "L4-A: bedtools genomecov"
if [ ! -f "${OUT}/coverage/genomecov.txt" ]; then
    bedtools genomecov -ibam "${BAM}" > "${OUT}/coverage/genomecov.txt"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4-B: deeptools bamCoverage (normalized signal)
# ===========================================================================
log_step "L4-B: deeptools bamCoverage"
if [ ! -f "${OUT}/coverage/sample.bw" ]; then
    bamCoverage --bam "${BAM}" --outFileName "${OUT}/coverage/sample.bw" \
                --binSize 100 --normalizeUsing RPGC \
                --effectiveGenomeSize 2900000 \
                --numberOfProcessors ${THREADS} 2>/dev/null || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5: MERGE
# ===========================================================================
log_step "L5-MERGE"

# Parse samtools stats
TOTAL_READS=$(grep "^SN" "${OUT}/stats/samtools_stats.txt" | grep "raw total sequences" | cut -f3)
MAPPED_READS=$(grep "^SN" "${OUT}/stats/samtools_stats.txt" | grep "reads mapped:" | cut -f3)
MAPPED_PCT=$(grep "mapped (" "${OUT}/stats/flagstat.txt" | head -1 | grep -oP '[\d.]+%' | tr -d '%')
PROPERLY_PAIRED=$(grep "properly paired" "${OUT}/stats/flagstat.txt" | grep -oP '[\d.]+%' | tr -d '%')
AVG_QUALITY=$(grep "^SN" "${OUT}/stats/samtools_stats.txt" | grep "average quality" | cut -f3)
AVG_INSERT=$(grep "^SN" "${OUT}/stats/samtools_stats.txt" | grep "insert size average" | cut -f3)
ERROR_RATE=$(grep "^SN" "${OUT}/stats/samtools_stats.txt" | grep "error rate" | cut -f3)

# Parse mosdepth
MEAN_COV=$(awk '$1=="total" {print $4}' "${OUT}/coverage/sample.mosdepth.summary.txt")
MIN_COV=$(awk '$1=="total" {print $5}' "${OUT}/coverage/sample.mosdepth.summary.txt")
MAX_COV=$(awk '$1=="total" {print $6}' "${OUT}/coverage/sample.mosdepth.summary.txt")

# Coverage breadth from genomecov (% bases with >=1x coverage)
BREADTH=$(awk '$1=="genome" && $2==0 {print 100-$5*100}' "${OUT}/coverage/genomecov.txt" | head -1)

cat > "${RES}/mapping_qc_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
mapped_reads,${MAPPED_READS}
mapping_rate,${MAPPED_PCT}
properly_paired_rate,${PROPERLY_PAIRED}
average_base_quality,${AVG_QUALITY}
average_insert_size,${AVG_INSERT}
error_rate,${ERROR_RATE}
mean_coverage,${MEAN_COV}
min_coverage,${MIN_COV}
max_coverage,${MAX_COV}
coverage_breadth_pct,${BREADTH}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/mapping_qc_report.csv"
echo ""
ls -lh "${RES}/"
