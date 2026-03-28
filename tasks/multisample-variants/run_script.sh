#!/bin/bash
set -e

# =============================================================================
# Task 24: Multi-sample Variant Calling and Comparison
#
# DAG (depth 6, per-sample parallel then merge):
#
# L0: reads_A + reads_B + reference
# L1: fastp A          fastp B
# L2: bwa A            bwa B
# L3: freebayes A      freebayes B
# L4: bcftools filter A  bcftools filter B
#     └──────────┬───────────┘
# L5: bcftools isec (shared vs unique variants)
#     └── variant annotation + summary
# L6: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

REFERENCE="${REF}/reference.fna"
SAMPLES=("sampleA" "sampleB")

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,aligned,variants,filtered,comparison} "${RES}"

# Per-sample processing (L1-L4)
for SAMPLE in "${SAMPLES[@]}"; do
    # L1: Trim
    log_step "L1: fastp ${SAMPLE}"
    if [ ! -f "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" ]; then
        fastp --in1 "${DATA}/${SAMPLE}_R1.fastq.gz" --in2 "${DATA}/${SAMPLE}_R2.fastq.gz" \
              --out1 "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" --out2 "${OUT}/trimmed/${SAMPLE}_R2.fastq.gz" \
              --detect_adapter_for_pe --thread ${THREADS} \
              --json "${OUT}/trimmed/${SAMPLE}_fastp.json"
    else echo "Skipping (exists)"; fi

    # L2: Align
    log_step "L2: bwa ${SAMPLE}"
    if [ ! -f "${OUT}/aligned/${SAMPLE}.bam" ]; then
        [ ! -f "${REFERENCE}.bwt" ] && bwa index "${REFERENCE}" 2>/dev/null
        bwa mem -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
                "${REFERENCE}" "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" "${OUT}/trimmed/${SAMPLE}_R2.fastq.gz" \
            | samtools sort -@ ${THREADS} -o "${OUT}/aligned/${SAMPLE}.bam"
        samtools index "${OUT}/aligned/${SAMPLE}.bam"
    else echo "Skipping (exists)"; fi

    # L3: Call variants
    log_step "L3: freebayes ${SAMPLE}"
    if [ ! -f "${OUT}/variants/${SAMPLE}.vcf" ]; then
        freebayes -f "${REFERENCE}" "${OUT}/aligned/${SAMPLE}.bam" \
                  --min-mapping-quality 20 --min-base-quality 20 \
                  > "${OUT}/variants/${SAMPLE}.vcf"
    else echo "Skipping (exists)"; fi

    # L4: Filter
    log_step "L4: bcftools filter ${SAMPLE}"
    if [ ! -f "${OUT}/filtered/${SAMPLE}.vcf.gz" ]; then
        bcftools filter -i 'QUAL>20 && INFO/DP>5' "${OUT}/variants/${SAMPLE}.vcf" \
            | bcftools view -Oz -o "${OUT}/filtered/${SAMPLE}.vcf.gz"
        bcftools index "${OUT}/filtered/${SAMPLE}.vcf.gz"
    else echo "Skipping (exists)"; fi
done

# L5: Compare variants between samples
log_step "L5: bcftools isec (shared vs unique)"
if [ ! -d "${OUT}/comparison/isec" ]; then
    bcftools isec -p "${OUT}/comparison/isec" \
                  "${OUT}/filtered/sampleA.vcf.gz" \
                  "${OUT}/filtered/sampleB.vcf.gz"
else echo "Skipping (exists)"; fi

# Parse comparison results
UNIQUE_A=$(grep -cv "^#" "${OUT}/comparison/isec/0000.vcf" 2>/dev/null || echo "0")
UNIQUE_B=$(grep -cv "^#" "${OUT}/comparison/isec/0001.vcf" 2>/dev/null || echo "0")
SHARED=$(grep -cv "^#" "${OUT}/comparison/isec/0002.vcf" 2>/dev/null || echo "0")

# Per-sample stats
VARS_A=$(bcftools stats "${OUT}/filtered/sampleA.vcf.gz" | grep "^SN" | grep "number of records" | cut -f4)
VARS_B=$(bcftools stats "${OUT}/filtered/sampleB.vcf.gz" | grep "^SN" | grep "number of records" | cut -f4)
SNPS_A=$(bcftools stats "${OUT}/filtered/sampleA.vcf.gz" | grep "^SN" | grep "number of SNPs" | cut -f4)
SNPS_B=$(bcftools stats "${OUT}/filtered/sampleB.vcf.gz" | grep "^SN" | grep "number of SNPs" | cut -f4)

# L6: MERGE
log_step "L6-MERGE"

cat > "${RES}/variant_comparison.csv" << CSVEOF
metric,value
total_variants_sampleA,${VARS_A}
total_variants_sampleB,${VARS_B}
snps_sampleA,${SNPS_A}
snps_sampleB,${SNPS_B}
unique_to_sampleA,${UNIQUE_A}
unique_to_sampleB,${UNIQUE_B}
shared_variants,${SHARED}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/variant_comparison.csv"
