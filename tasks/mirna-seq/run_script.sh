#!/bin/bash
set -euo pipefail

# =============================================================================
# Task: Small RNA-seq / miRNA Discovery
#
# DAG structure (depth 7, 2 convergence points):
#
# L0: SE reads (small RNA-seq)
# L1: fastp (adapter trim — critical for 18-25nt miRNAs)
# L2: bowtie (align to mature miRNA reference)
#     ├─────────────────────────────────────────┐
# L3: samtools (count mature hits)        bowtie (align to hairpin reference)
#     │                                         │
# L4: mirtop (standardize miRNA counts) ◄──────┘  [CONVERGENCE 1: merge mature+hairpin]
#     │                      │
# L5: mirtrace (QC + origin) │
#     │                      │
# L6: count summary    isomiR analysis
#     └──────────┬───────────┘
# L7: MERGE ◄────────────────────────────────────  [CONVERGENCE 2: QC + counts]
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

MATURE="${REF}/mature_human.fa"
HAIRPIN="${REF}/hairpin_human.fa"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{trimmed,mature_align,hairpin_align,mirtop,mirtrace,counts} "${RES}"

# L1: Adapter trimming (critical — adapters dominate miRNA reads)
log_step "L1: fastp adapter trim"
if [ ! -f "${OUT}/trimmed/reads.fastq.gz" ]; then
    # Trim TruSeq Small RNA 3' adapter
    fastp --in1 "${DATA}/reads.fastq.gz" \
          --out1 "${OUT}/trimmed/reads.fastq.gz" \
          --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
          --length_required 18 --length_limit 30 \
          --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
    zcat "${OUT}/trimmed/reads.fastq.gz" > "${OUT}/trimmed/reads.fastq"
fi

# L2: Build bowtie indexes + align to mature miRNAs
log_step "L2: bowtie index + align to mature"
if [ ! -f "${OUT}/mature_align/aligned.bam" ]; then
    # Replace U with T in miRBase FASTA (RNA→DNA)
    sed 's/U/T/g' "${MATURE}" > "${OUT}/mature_align/mature_dna.fa"
    bowtie-build "${OUT}/mature_align/mature_dna.fa" "${OUT}/mature_align/mature_idx" --quiet
    bowtie -x "${OUT}/mature_align/mature_idx" \
           -q "${OUT}/trimmed/reads.fastq" \
           -v 1 -k 5 --best --strata -p ${THREADS} -S 2>/dev/null \
        | samtools view -bS -F 4 | samtools sort -o "${OUT}/mature_align/aligned.bam"
    samtools index "${OUT}/mature_align/aligned.bam"
fi

# L3 RIGHT: Align to hairpin
log_step "L3: bowtie align to hairpin"
if [ ! -f "${OUT}/hairpin_align/aligned.bam" ]; then
    sed 's/U/T/g' "${HAIRPIN}" > "${OUT}/hairpin_align/hairpin_dna.fa"
    bowtie-build "${OUT}/hairpin_align/hairpin_dna.fa" "${OUT}/hairpin_align/hairpin_idx" --quiet
    bowtie -x "${OUT}/hairpin_align/hairpin_idx" \
           -q "${OUT}/trimmed/reads.fastq" \
           -v 1 -k 5 --best --strata -p ${THREADS} -S 2>/dev/null \
        | samtools view -bS -F 4 | samtools sort -o "${OUT}/hairpin_align/aligned.bam"
    samtools index "${OUT}/hairpin_align/aligned.bam"
fi

# L3 LEFT: Count mature hits
log_step "L3: count mature miRNA hits"
samtools idxstats "${OUT}/mature_align/aligned.bam" | \
    awk '$3 > 0 {print $1"\t"$3}' | sort -k2 -rn > "${OUT}/counts/mature_counts.tsv"

# L4: mirtrace QC
log_step "L4: mirtrace QC"
if [ ! -d "${OUT}/mirtrace/mirtrace-results" ]; then
    mirtrace qc --species hsa --protocol illumina \
             --output-dir "${OUT}/mirtrace" \
             "${OUT}/trimmed/reads.fastq.gz" 2>&1 || true
fi

# L5-L7: MERGE
log_step "MERGE"

# Parse fastp stats
TOTAL_READS=$(python3 -c "import json; d=json.load(open('${OUT}/trimmed/fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
PASSED_READS=$(python3 -c "import json; d=json.load(open('${OUT}/trimmed/fastp.json')); print(d['summary']['after_filtering']['total_reads'])")

# Alignment stats
MATURE_MAPPED=$(samtools view -c -F 4 "${OUT}/mature_align/aligned.bam" 2>/dev/null || echo "0")
HAIRPIN_MAPPED=$(samtools view -c -F 4 "${OUT}/hairpin_align/aligned.bam" 2>/dev/null || echo "0")

# Unique miRNAs detected
UNIQUE_MIRNAS=$(wc -l < "${OUT}/counts/mature_counts.tsv" | tr -d ' ')

# Top miRNA
TOP_MIRNA=$(head -1 "${OUT}/counts/mature_counts.tsv" | cut -f1 || echo "N/A")
TOP_COUNT=$(head -1 "${OUT}/counts/mature_counts.tsv" | cut -f2 || echo "0")

# mirtrace QC
MIRNA_PCT=$(grep "miRNA" "${OUT}/mirtrace/mirtrace-results"/*.html 2>/dev/null | grep -oP '[\d.]+%' | head -1 | tr -d '%' || echo "N/A")

cat > "${RES}/mirna_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
passed_filter,${PASSED_READS}
mature_mapped_reads,${MATURE_MAPPED}
hairpin_mapped_reads,${HAIRPIN_MAPPED}
unique_mirnas_detected,${UNIQUE_MIRNAS}
top_mirna,${TOP_MIRNA}
top_mirna_count,${TOP_COUNT}
CSVEOF

# Copy count table
cp "${OUT}/counts/mature_counts.tsv" "${RES}/"

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/mirna_report.csv"
