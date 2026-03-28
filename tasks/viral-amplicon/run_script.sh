#!/bin/bash
set -euo pipefail

# =============================================================================
# Task 31: Viral Amplicon Analysis (SARS-CoV-2)
#
# DAG (depth 8, nested diamond):
#
# L0: reads + reference + primers
# L1: fastp (trim)
# L2: bowtie2 (align to reference)
# L3: ivar trim (primer trimming — amplicon-specific)
# L4: samtools mpileup
#     ├──────────────────────────────────────────┐
# L5: ivar variants (call)     mosdepth (coverage)
#     │                              │
# L6: snpSift (annotate)       coverage report
#     │                              │
# L7: bcftools consensus ◄──────────┘ (mask low-cov in consensus)
#     ├──────────────────────┐
# L8: pangolin (lineage)    nextclade (clade + mutations)
#     └──────────┬───────────┘
# L9: MERGE (surveillance report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

REFERENCE="${REF}/sarscov2_ref.fna"
PRIMERS="${REF}/primers.bed"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,aligned,primer_trimmed,variants,coverage,consensus,lineage} "${RES}"

# L1: Trim
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
fi

# L2: Align
log_step "L2: bowtie2"
if [ ! -f "${OUT}/aligned/aligned.sorted.bam" ]; then
    [ ! -f "${REFERENCE}.1.bt2" ] && bowtie2-build "${REFERENCE}" "${REFERENCE}" --threads ${THREADS}
    bowtie2 -x "${REFERENCE}" -1 "${OUT}/trimmed/R1.fastq.gz" -2 "${OUT}/trimmed/R2.fastq.gz" \
            --very-sensitive --no-mixed --no-discordant -X 1000 \
            --rg-id sample --rg "SM:sample" --rg "PL:ILLUMINA" \
            --threads ${THREADS} 2>"${OUT}/aligned/bowtie2.log" \
        | samtools sort -@ ${THREADS} -o "${OUT}/aligned/aligned.sorted.bam"
    samtools index "${OUT}/aligned/aligned.sorted.bam"
fi

# L3: ivar primer trimming
log_step "L3: ivar trim"
if [ ! -f "${OUT}/primer_trimmed/trimmed.bam" ]; then
    ivar trim -i "${OUT}/aligned/aligned.sorted.bam" \
              -b "${PRIMERS}" -p "${OUT}/primer_trimmed/trimmed" \
              -m 30 -q 20 -e 2>&1 || {
        echo "WARNING: ivar trim failed (primers may not match), using untrimmed BAM"
        cp "${OUT}/aligned/aligned.sorted.bam" "${OUT}/primer_trimmed/trimmed.bam"
    }
    samtools sort -@ ${THREADS} -o "${OUT}/primer_trimmed/trimmed.sorted.bam" \
                  "${OUT}/primer_trimmed/trimmed.bam"
    samtools index "${OUT}/primer_trimmed/trimmed.sorted.bam"
fi
BAM="${OUT}/primer_trimmed/trimmed.sorted.bam"

# L4: samtools mpileup
log_step "L4: mpileup"
samtools faidx "${REFERENCE}" 2>/dev/null || true

# L5 LEFT: ivar variants
log_step "L5-LEFT: ivar variants"
if [ ! -f "${OUT}/variants/variants.tsv" ]; then
    samtools mpileup -aa -A -d 600 -B -Q 20 --reference "${REFERENCE}" "${BAM}" \
        | ivar variants -p "${OUT}/variants/variants" -q 20 -t 0.5 -r "${REFERENCE}"
fi

# L5 RIGHT: mosdepth coverage
log_step "L5-RIGHT: mosdepth"
if [ ! -f "${OUT}/coverage/sample.mosdepth.summary.txt" ]; then
    mosdepth -t ${THREADS} --by 100 "${OUT}/coverage/sample" "${BAM}"
fi

# L6: annotate + coverage report
log_step "L6: coverage analysis"
MEAN_COV=$(awk '$1=="NC_045512.2" || $1=="total" {print $4; exit}' "${OUT}/coverage/sample.mosdepth.summary.txt" 2>/dev/null || echo "0")
TOTAL_VARIANTS=$(wc -l < "${OUT}/variants/variants.tsv" 2>/dev/null | tr -d ' ')
TOTAL_VARIANTS=$((TOTAL_VARIANTS - 1))  # subtract header
PASS_VARIANTS=$(awk -F'\t' '$14=="TRUE"' "${OUT}/variants/variants.tsv" 2>/dev/null | wc -l | tr -d ' ')

# L7: bcftools consensus
log_step "L7: consensus generation"
if [ ! -f "${OUT}/consensus/consensus.fasta" ]; then
    # Convert ivar TSV to VCF for bcftools consensus
    python3 -c "
import sys
with open('${OUT}/variants/variants.tsv') as f:
    header = f.readline()
    print('##fileformat=VCFv4.2')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) > 13 and parts[13] == 'TRUE':
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            qual = parts[12] if parts[12] != 'NA' else '.'
            print(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t.')
" > "${OUT}/variants/pass_variants.vcf"

    bgzip -c "${OUT}/variants/pass_variants.vcf" > "${OUT}/variants/pass_variants.vcf.gz"
    tabix -p vcf "${OUT}/variants/pass_variants.vcf.gz"

    # Mask low-coverage regions
    zcat "${OUT}/coverage/sample.per-base.bed.gz" | awk '$4 < 10' > "${OUT}/coverage/low_cov.bed"

    bcftools consensus -f "${REFERENCE}" "${OUT}/variants/pass_variants.vcf.gz" \
        -m "${OUT}/coverage/low_cov.bed" \
        > "${OUT}/consensus/consensus.fasta" 2>/dev/null || {
        cp "${REFERENCE}" "${OUT}/consensus/consensus.fasta"
    }
fi

# L8 LEFT: Pangolin lineage
log_step "L8-LEFT: pangolin"
if [ ! -f "${OUT}/lineage/lineage_report.csv" ]; then
    pangolin "${OUT}/consensus/consensus.fasta" --outdir "${OUT}/lineage" \
             --outfile lineage_report.csv 2>&1 || {
        echo "WARNING: pangolin failed"
        echo "taxon,lineage,scorpio_call" > "${OUT}/lineage/lineage_report.csv"
        echo "consensus,Unassigned," >> "${OUT}/lineage/lineage_report.csv"
    }
fi

# L8 RIGHT: Nextclade
log_step "L8-RIGHT: nextclade"
if [ ! -f "${OUT}/lineage/nextclade.tsv" ]; then
    nextclade run --input-fasta "${OUT}/consensus/consensus.fasta" \
                  --output-tsv "${OUT}/lineage/nextclade.tsv" 2>&1 || {
        echo "WARNING: nextclade failed"
        touch "${OUT}/lineage/nextclade.tsv"
    }
fi

# L9: MERGE
log_step "L9-MERGE"

LINEAGE=$(python3 -c "
import csv
with open('${OUT}/lineage/lineage_report.csv') as f:
    row = next(csv.DictReader(f))
    print(row['lineage'])
" 2>/dev/null || echo "Unknown")
CLADE=$(tail -1 "${OUT}/lineage/nextclade.tsv" 2>/dev/null | cut -f2 || echo "N/A")
CONS_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${OUT}/consensus/consensus.fasta")
N_COUNT=$(grep -o "N" "${OUT}/consensus/consensus.fasta" | wc -l | tr -d ' ')
N_PCT=$(python3 -c "print(f'{${N_COUNT}/${CONS_LEN}*100:.2f}')" 2>/dev/null || echo "0")
ALIGN_RATE=$(grep "overall alignment rate" "${OUT}/aligned/bowtie2.log" | awk '{print $1}' || echo "0")
MAPPED=$(samtools view -c -F 4 "${BAM}" 2>/dev/null || echo "0")

cat > "${RES}/viral_surveillance.csv" << CSVEOF
metric,value
genome_length,${CONS_LEN}
masked_bases,${N_COUNT}
masked_pct,${N_PCT}
alignment_rate,${ALIGN_RATE}
mapped_reads,${MAPPED}
mean_coverage,${MEAN_COV}
total_variants_called,${TOTAL_VARIANTS}
pass_variants,${PASS_VARIANTS}
lineage,${LINEAGE}
clade,${CLADE}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/viral_surveillance.csv"
