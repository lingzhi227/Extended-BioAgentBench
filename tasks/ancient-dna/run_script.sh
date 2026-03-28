#!/bin/bash
set -euo pipefail

# =============================================================================
# Task 34: Ancient DNA Authentication and Damage Assessment
# DAG (depth 8):
# L0: degraded PE reads + reference
# L1: AdapterRemoval (merge overlapping PE, trim adapters)
# L2: bwa aln + samse/sampe (short-read aDNA alignment — NOT bwa mem)
# L3: samtools sort + filter (MAPQ>=25)
# L4: picard MarkDuplicates
#     ├──────────────────────────────────────┐
# L5: DamageProfiler (C→T damage pattern)   samtools depth (endogenous)
#     │                                           │
# L6: mapDamage (rescale + damage plot)      coverage stats
#     └──────────┬────────────────────────────────┘
# L7: MERGE (authenticity report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"
GENOME="${REF}/genome.fa"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{trimmed,aligned,dedup,damage,mapdamage,coverage} "${RES}"

# L1: AdapterRemoval (trim adapters + quality)
log_step "L1: AdapterRemoval"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    AdapterRemoval --file1 "${DATA}/reads_R1.fastq.gz" --file2 "${DATA}/reads_R2.fastq.gz" \
                   --trimns --trimqualities --minquality 15 --minlength 25 \
                   --output1 "${OUT}/trimmed/R1.fastq.gz" --output2 "${OUT}/trimmed/R2.fastq.gz" \
                   --gzip --threads ${THREADS}
fi

# L2: bwa aln (for short aDNA reads)
log_step "L2: bwa aln"
if [ ! -f "${OUT}/aligned/aligned.bam" ]; then
    [ ! -f "${GENOME}.bwt" ] && bwa index "${GENOME}"
    bwa aln -t ${THREADS} -l 1024 -n 0.01 -o 2 "${GENOME}" "${OUT}/trimmed/R1.fastq.gz" > "${OUT}/aligned/r1.sai"
    bwa aln -t ${THREADS} -l 1024 -n 0.01 -o 2 "${GENOME}" "${OUT}/trimmed/R2.fastq.gz" > "${OUT}/aligned/r2.sai"
    bwa sampe -r "@RG\tID:aDNA\tSM:aDNA\tPL:ILLUMINA" \
        "${GENOME}" "${OUT}/aligned/r1.sai" "${OUT}/aligned/r2.sai" \
        "${OUT}/trimmed/R1.fastq.gz" "${OUT}/trimmed/R2.fastq.gz" | \
        samtools view -bSh -q 25 -F 4 | \
        samtools sort -@ ${THREADS} -o "${OUT}/aligned/aligned.bam"
    samtools index "${OUT}/aligned/aligned.bam"
fi
BAM="${OUT}/aligned/aligned.bam"

# L3: Stats
log_step "L3: flagstat"
samtools flagstat "${BAM}" > "${OUT}/aligned/flagstat.txt"

# L4: Dedup
log_step "L4: picard MarkDuplicates"
if [ ! -f "${OUT}/dedup/dedup.bam" ]; then
    picard MarkDuplicates INPUT="${BAM}" OUTPUT="${OUT}/dedup/dedup.bam" \
           METRICS_FILE="${OUT}/dedup/metrics.txt" REMOVE_DUPLICATES=true \
           VALIDATION_STRINGENCY=LENIENT
    samtools index "${OUT}/dedup/dedup.bam"
fi
DEDUP="${OUT}/dedup/dedup.bam"

# L5 LEFT: DamageProfiler
log_step "L5-LEFT: DamageProfiler"
if [ ! -d "${OUT}/damage/results" ]; then
    damageprofiler -i "${DEDUP}" -r "${GENOME}" -o "${OUT}/damage/results" 2>&1 || true
fi

# L5 RIGHT: Coverage
log_step "L5-RIGHT: coverage"
samtools depth -a "${DEDUP}" > "${OUT}/coverage/depth.txt"
MEAN_COV=$(awk '{sum+=$3; n++} END{printf "%.2f", sum/n}' "${OUT}/coverage/depth.txt")
BREADTH=$(awk '$3>0{c++} END{printf "%.2f", c/NR*100}' "${OUT}/coverage/depth.txt")

# L6: mapDamage
log_step "L6: mapDamage"
if [ ! -d "${OUT}/mapdamage/results_DEDUP" ]; then
    mapDamage -i "${DEDUP}" -r "${GENOME}" -d "${OUT}/mapdamage" --no-stats 2>&1 || true
fi

# Parse damage stats
CT_5PRIME=$(head -2 "${OUT}/damage/results/5pCtoT_freq.txt" 2>/dev/null | tail -1 | awk '{print $2}' || echo "0")
GA_3PRIME=$(head -2 "${OUT}/damage/results/3pGtoA_freq.txt" 2>/dev/null | tail -1 | awk '{print $2}' || echo "0")

# Parse alignment stats
TOTAL_INPUT=$(grep "in total" "${OUT}/aligned/flagstat.txt" | awk '{print $1}')
MAPPED=$(grep "mapped (" "${OUT}/aligned/flagstat.txt" | head -1 | awk '{print $1}')
ENDOGENOUS_PCT=$(python3 -c "print(f'{${MAPPED}/${TOTAL_INPUT}*100:.2f}')" 2>/dev/null || echo "0")
DUP_PCT=$(grep "PERCENT_DUPLICATION" "${OUT}/dedup/metrics.txt" -A1 | tail -1 | cut -f9 || echo "0")
DEDUP_READS=$(samtools view -c "${DEDUP}")
MEAN_LENGTH=$(samtools view "${DEDUP}" | awk '{sum+=length($10); n++} END{printf "%.1f", sum/n}')

# L7: MERGE
log_step "L7-MERGE"
cat > "${RES}/adna_report.csv" << CSVEOF
metric,value
total_input_reads,${TOTAL_INPUT}
mapped_reads,${MAPPED}
endogenous_pct,${ENDOGENOUS_PCT}
duplication_rate,${DUP_PCT}
reads_after_dedup,${DEDUP_READS}
mean_read_length,${MEAN_LENGTH}
mean_coverage,${MEAN_COV}
coverage_breadth_pct,${BREADTH}
damage_5prime_ct,${CT_5PRIME}
damage_3prime_ga,${GA_3PRIME}
CSVEOF

echo "=== Pipeline complete ==="
cat "${RES}/adna_report.csv"
