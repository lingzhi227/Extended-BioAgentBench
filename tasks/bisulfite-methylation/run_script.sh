#!/bin/bash
set -euo pipefail

# =============================================================================
# Task 32: Bisulfite Methylation Analysis
#
# DAG (depth 8):
# L0: WGBS reads + genome
# L1: Trim Galore (adapter + quality trim)
# L2: bismark_genome_preparation (index)
# L3: bismark align (bisulfite-aware alignment)
# L4: deduplicate_bismark (remove PCR dups)
# L5: bismark_methylation_extractor (extract CpG methylation)
#     ├──────────────────────────────────────────┐
# L6: ├── bedGraph (per-CpG methylation levels)  │
#     └── M-bias plot (strand bias QC)           │
#          │                                      │
# L7: ├── coverage2cytosine (genome-wide report) │
#     └── bismark2summary (comprehensive QC)     │
# L8: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME_DIR="${REF}"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,aligned,dedup,methylation,reports} "${RES}"

# L1: Trim
log_step "L1: Trim Galore"
if [ ! -f "${OUT}/trimmed/reads_R1_val_1.fq.gz" ]; then
    trim_galore --paired --fastqc -o "${OUT}/trimmed" \
                "${DATA}/reads_R1.fastq.gz" "${DATA}/reads_R2.fastq.gz"
fi

# L2: Prepare bisulfite genome
log_step "L2: bismark genome preparation"
if [ ! -d "${GENOME_DIR}/Bisulfite_Genome" ]; then
    bismark_genome_preparation "${GENOME_DIR}"
fi

# L3: Bismark alignment
log_step "L3: bismark align"
if [ ! -f "${OUT}/aligned/reads_R1_val_1_bismark_bt2_pe.bam" ]; then
    bismark --genome "${GENOME_DIR}" \
            -1 "${OUT}/trimmed/reads_R1_val_1.fq.gz" \
            -2 "${OUT}/trimmed/reads_R2_val_2.fq.gz" \
            --output_dir "${OUT}/aligned" \
            --parallel $((THREADS / 2 > 0 ? THREADS / 2 : 1))
fi
BAM=$(ls "${OUT}/aligned/"*_bismark_bt2_pe.bam 2>/dev/null | head -1)

# L4: Dedup
log_step "L4: deduplicate_bismark"
if [ ! -f "${OUT}/dedup/$(basename ${BAM%.bam}).deduplicated.bam" ]; then
    deduplicate_bismark --bam "${BAM}" --output_dir "${OUT}/dedup"
fi
DEDUP_BAM=$(ls "${OUT}/dedup/"*.deduplicated.bam 2>/dev/null | head -1)

# L5: Extract methylation
log_step "L5: bismark_methylation_extractor"
if [ ! -f "${OUT}/methylation/done" ]; then
    bismark_methylation_extractor --paired-end --bedGraph --cytosine_report \
        --genome_folder "${GENOME_DIR}" \
        --output "${OUT}/methylation" \
        "${DEDUP_BAM}"
    touch "${OUT}/methylation/done"
fi

# L6: Reports
log_step "L6: bismark reports"
bismark2report --alignment_report "${OUT}/aligned/"*PE_report.txt \
               --dedup_report "${OUT}/dedup/"*deduplication_report.txt \
               --mbias_report "${OUT}/methylation/"*M-bias.txt \
               --dir "${OUT}/reports" 2>/dev/null || true

bismark2summary "${OUT}/aligned/"*PE_report.txt 2>/dev/null || true

# L7-L8: MERGE
log_step "MERGE"

# Parse alignment report
ALIGN_REPORT=$(ls ${OUT}/aligned/*PE_report.txt 2>/dev/null | head -1)
TOTAL_PAIRS=$(grep "Sequence pairs analysed" "${ALIGN_REPORT}" 2>/dev/null | awk -F'\t' '{print $NF}' | tr -d ' ' || echo "0")
ALIGNED_PAIRS=$(grep "Number of paired-end alignments" "${ALIGN_REPORT}" 2>/dev/null | awk -F'\t' '{print $NF}' | tr -d ' ' || echo "0")
MAPPING_EFF=$(grep "Mapping efficiency" "${ALIGN_REPORT}" 2>/dev/null | awk -F: '{print $2}' | tr -d ' \t%' || echo "0")

# Parse dedup report
DEDUP_REPORT=$(ls ${OUT}/dedup/*deduplication_report.txt 2>/dev/null | head -1)
DUP_PCT=$(grep "Total duplicate" "${DEDUP_REPORT}" 2>/dev/null | grep -oP '[\d.]+%' | tr -d '%' || echo "0")

# Parse methylation
METHYL_CpG=$(grep "C methylated in CpG context" "${ALIGN_REPORT}" 2>/dev/null | awk -F: '{print $2}' | tr -d ' \t%' || echo "0")
METHYL_CHG=$(grep "C methylated in CHG context" "${ALIGN_REPORT}" 2>/dev/null | awk -F: '{print $2}' | tr -d ' \t%' || echo "0")
METHYL_CHH=$(grep "C methylated in CHH context" "${ALIGN_REPORT}" 2>/dev/null | awk -F: '{print $2}' | tr -d ' \t%' || echo "0")

# Count CpG sites from cytosine report
CYT_REPORT=$(ls "${OUT}/methylation/"*.CpG_report.txt 2>/dev/null | head -1)
if [ -n "$CYT_REPORT" ] && [ -f "$CYT_REPORT" ]; then
    TOTAL_CPG=$(wc -l < "$CYT_REPORT" | tr -d ' ')
    METHYLATED_CPG=$(awk '$4 > 0' "$CYT_REPORT" | wc -l | tr -d ' ')
else
    TOTAL_CPG=0; METHYLATED_CPG=0
fi

cat > "${RES}/methylation_report.csv" << CSVEOF
metric,value
total_read_pairs,${TOTAL_PAIRS}
aligned_pairs,${ALIGNED_PAIRS}
mapping_efficiency,${MAPPING_EFF}
duplication_rate,${DUP_PCT}
cpg_methylation_pct,${METHYL_CpG}
chg_methylation_pct,${METHYL_CHG}
chh_methylation_pct,${METHYL_CHH}
total_cpg_sites,${TOTAL_CPG}
methylated_cpg_sites,${METHYLATED_CPG}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/methylation_report.csv"
