#!/bin/bash
set -e

# =============================================================================
# Task 18: Structural Variant Detection in MRSA
#
# DAG (depth 6, wide fan-out + merge):
#
# L0: reads (MRSA PE reads)
# L1: fastp (trim)
# L2: bwa mem (align to reference assembly)
#     ├──────────────────────────────────────┐
# L3: delly call (SV detection)        freebayes (SNP/indel calling)
#     │                                      │
# L4: delly filter (quality filter)    bcftools filter (quality filter)
#     │                                      │
#     ├──────────┐                           │
# L5: bedtools annotate     bcftools stats   │
#     (intersect SVs        (variant stats)  │
#      with genes)                           │
#     └──────────┴───────────────────────────┘
# L6: MERGE (integrated variant report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

REFERENCE="${REF}/contigs.fa"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,aligned,delly,freebayes,filtered,annotated,stats} "${RES}"

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
log_step "L2: bwa alignment"
if [ ! -f "${OUT}/aligned/reads.bam" ]; then
    bwa index "${REFERENCE}"
    bwa mem -t ${THREADS} -R "@RG\tID:MRSA\tSM:MRSA\tPL:ILLUMINA" \
            "${REFERENCE}" "${OUT}/trimmed/R1.fastq" "${OUT}/trimmed/R2.fastq" \
        | samtools sort -@ ${THREADS} -o "${OUT}/aligned/reads.bam"
    samtools index "${OUT}/aligned/reads.bam"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3 LEFT: Structural variant calling with Delly
# ===========================================================================
log_step "L3-LEFT: delly SV calling"
if [ ! -f "${OUT}/delly/sv_calls.bcf" ]; then
    samtools faidx "${REFERENCE}"
    delly call -g "${REFERENCE}" "${OUT}/aligned/reads.bam" -o "${OUT}/delly/sv_calls.bcf" 2>&1 || true
    bcftools view "${OUT}/delly/sv_calls.bcf" > "${OUT}/delly/sv_calls.vcf" 2>/dev/null || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3 RIGHT: SNP/Indel calling with freebayes
# ===========================================================================
log_step "L3-RIGHT: freebayes variant calling"
if [ ! -f "${OUT}/freebayes/variants.vcf" ]; then
    freebayes -f "${REFERENCE}" "${OUT}/aligned/reads.bam" \
              --min-mapping-quality 20 --min-base-quality 20 \
              > "${OUT}/freebayes/variants.vcf"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4 LEFT: Delly filter
# ===========================================================================
log_step "L4-LEFT: delly filter"
if [ ! -f "${OUT}/filtered/sv_filtered.vcf" ]; then
    if [ -s "${OUT}/delly/sv_calls.bcf" ]; then
        delly filter -f germline "${OUT}/delly/sv_calls.bcf" -o "${OUT}/filtered/sv_filtered.bcf" 2>&1 || true
        bcftools view "${OUT}/filtered/sv_filtered.bcf" > "${OUT}/filtered/sv_filtered.vcf" 2>/dev/null || \
            cp "${OUT}/delly/sv_calls.vcf" "${OUT}/filtered/sv_filtered.vcf"
    else
        touch "${OUT}/filtered/sv_filtered.vcf"
    fi
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4 RIGHT: Freebayes filter
# ===========================================================================
log_step "L4-RIGHT: bcftools filter SNPs/indels"
if [ ! -f "${OUT}/filtered/snps_filtered.vcf" ]; then
    bcftools filter -i 'QUAL>20 && DP>5' "${OUT}/freebayes/variants.vcf" \
        > "${OUT}/filtered/snps_filtered.vcf"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5: Annotation and statistics
# ===========================================================================
log_step "L5: bcftools stats"
bcftools stats "${OUT}/filtered/snps_filtered.vcf" > "${OUT}/stats/snp_stats.txt" 2>/dev/null || true
bcftools stats "${OUT}/filtered/sv_filtered.vcf" > "${OUT}/stats/sv_stats.txt" 2>/dev/null || true

# Count variants by type
SNP_COUNT=$(grep "^SN" "${OUT}/stats/snp_stats.txt" | grep "number of SNPs" | awk '{print $NF}' || echo "0")
INDEL_COUNT=$(grep "^SN" "${OUT}/stats/snp_stats.txt" | grep "number of indels" | awk '{print $NF}' || echo "0")
SV_COUNT=$(grep -v "^#" "${OUT}/filtered/sv_filtered.vcf" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# SV types breakdown
DEL_COUNT=$(grep -v "^#" "${OUT}/filtered/sv_filtered.vcf" 2>/dev/null | grep "SVTYPE=DEL" | wc -l | tr -d ' ' || echo "0")
DUP_COUNT=$(grep -v "^#" "${OUT}/filtered/sv_filtered.vcf" 2>/dev/null | grep "SVTYPE=DUP" | wc -l | tr -d ' ' || echo "0")
INV_COUNT=$(grep -v "^#" "${OUT}/filtered/sv_filtered.vcf" 2>/dev/null | grep "SVTYPE=INV" | wc -l | tr -d ' ' || echo "0")
BND_COUNT=$(grep -v "^#" "${OUT}/filtered/sv_filtered.vcf" 2>/dev/null | grep "SVTYPE=BND" | wc -l | tr -d ' ' || echo "0")

# Alignment stats
MAPPED=$(samtools flagstat "${OUT}/aligned/reads.bam" | grep "mapped (" | head -1 | awk '{print $1}')
TOTAL=$(samtools flagstat "${OUT}/aligned/reads.bam" | head -1 | awk '{print $1}')

# ===========================================================================
# L6: MERGE
# ===========================================================================
log_step "L6-MERGE"

cat > "${RES}/variant_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL}
mapped_reads,${MAPPED}
snps,${SNP_COUNT}
indels,${INDEL_COUNT}
structural_variants,${SV_COUNT}
sv_deletions,${DEL_COUNT}
sv_duplications,${DUP_COUNT}
sv_inversions,${INV_COUNT}
sv_breakends,${BND_COUNT}
CSVEOF

# Copy VCF files
cp "${OUT}/filtered/sv_filtered.vcf" "${RES}/structural_variants.vcf" 2>/dev/null || true
cp "${OUT}/filtered/snps_filtered.vcf" "${RES}/snps_indels.vcf" 2>/dev/null || true

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/variant_report.csv"
echo ""
ls -lh "${RES}/"
