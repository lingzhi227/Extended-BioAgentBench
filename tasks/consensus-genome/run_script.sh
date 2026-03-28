#!/bin/bash
set -e

# =============================================================================
# Task 25: Consensus Genome Generation
#
# DAG (depth 6):
#
# L0: reads + reference
# L1: fastp (trim)
# L2: bwa mem (align)
# L3: ├── samtools mpileup → bcftools call (variant calling)
#     └── mosdepth (coverage for masking)
# L4: bcftools consensus (apply variants to reference)
#     + bedtools (mask low-coverage regions)
# L5: ├── seqkit stats (consensus QC)
#     └── nucmer (compare consensus vs reference)
# L6: MERGE
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

mkdir -p "${OUT}"/{trimmed,aligned,variants,coverage,consensus,comparison} "${RES}"

# L1
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
else echo "Skipping (exists)"; fi

# L2
log_step "L2: bwa"
if [ ! -f "${OUT}/aligned/reads.bam" ]; then
    [ ! -f "${REFERENCE}.bwt" ] && bwa index "${REFERENCE}" 2>/dev/null
    bwa mem -t ${THREADS} -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
            "${REFERENCE}" "${OUT}/trimmed/R1.fastq.gz" "${OUT}/trimmed/R2.fastq.gz" \
        | samtools sort -@ ${THREADS} -o "${OUT}/aligned/reads.bam"
    samtools index "${OUT}/aligned/reads.bam"
else echo "Skipping (exists)"; fi

# L3-A: Variant calling
log_step "L3-A: bcftools mpileup + call"
if [ ! -f "${OUT}/variants/calls.vcf.gz" ]; then
    samtools faidx "${REFERENCE}" 2>/dev/null || true
    bcftools mpileup -f "${REFERENCE}" "${OUT}/aligned/reads.bam" --max-depth 500 -Q 20 \
        | bcftools call -c --ploidy 1 -Oz -o "${OUT}/variants/calls.vcf.gz"
    bcftools index "${OUT}/variants/calls.vcf.gz"
else echo "Skipping (exists)"; fi

# L3-B: Coverage for masking
log_step "L3-B: mosdepth coverage"
if [ ! -f "${OUT}/coverage/sample.per-base.bed.gz" ]; then
    mosdepth -t ${THREADS} "${OUT}/coverage/sample" "${OUT}/aligned/reads.bam"
else echo "Skipping (exists)"; fi

# L4: Generate consensus with low-coverage masking
log_step "L4: consensus generation"
if [ ! -f "${OUT}/consensus/consensus.fasta" ]; then
    # Create mask for regions with <5x coverage
    zcat "${OUT}/coverage/sample.per-base.bed.gz" | awk '$4 < 5' > "${OUT}/coverage/low_cov.bed"

    # Apply variants and mask
    bcftools consensus -f "${REFERENCE}" "${OUT}/variants/calls.vcf.gz" \
        -m "${OUT}/coverage/low_cov.bed" \
        > "${OUT}/consensus/consensus.fasta"
else echo "Skipping (exists)"; fi

# L5-A: Consensus QC
log_step "L5-A: seqkit stats"
seqkit stats "${OUT}/consensus/consensus.fasta" -T > "${OUT}/consensus/stats.tsv"

# L5-B: Compare consensus vs reference
log_step "L5-B: nucmer comparison"
if [ ! -f "${OUT}/comparison/dnadiff.report" ]; then
    nucmer --prefix="${OUT}/comparison/alignment" "${REFERENCE}" "${OUT}/consensus/consensus.fasta"
    dnadiff -d "${OUT}/comparison/alignment.delta" -p "${OUT}/comparison/dnadiff"
else echo "Skipping (exists)"; fi

# L6: MERGE
log_step "L6-MERGE"

CONS_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${OUT}/consensus/consensus.fasta")
REF_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${REFERENCE}")
N_COUNT=$(grep -o "N" "${OUT}/consensus/consensus.fasta" | wc -l | tr -d ' ')
N_PCT=$(python3 -c "print(f'{${N_COUNT}/${CONS_LEN}*100:.2f}')")
IDENTITY=$(grep "AvgIdentity" "${OUT}/comparison/dnadiff.report" | head -1 | awk '{print $2}')
ALIGNED=$(grep "AlignedBases" "${OUT}/comparison/dnadiff.report" | head -1 | awk '{print $2}' | sed 's/(.*//')
SNPS_VS_REF=$(grep "TotalSNPs" "${OUT}/comparison/dnadiff.report" | head -1 | awk '{print $2}')
MEAN_COV=$(awk '$1=="total" {print $4}' "${OUT}/coverage/sample.mosdepth.summary.txt")
LOW_COV_BASES=$(wc -l < "${OUT}/coverage/low_cov.bed" | tr -d ' ')

cat > "${RES}/consensus_report.csv" << CSVEOF
metric,value
reference_length,${REF_LEN}
consensus_length,${CONS_LEN}
masked_bases_n,${N_COUNT}
masked_pct,${N_PCT}
mean_coverage,${MEAN_COV}
low_coverage_regions,${LOW_COV_BASES}
identity_to_reference,${IDENTITY}
aligned_bases,${ALIGNED}
snps_vs_reference,${SNPS_VS_REF}
CSVEOF

cp "${OUT}/consensus/consensus.fasta" "${RES}/"

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/consensus_report.csv"
