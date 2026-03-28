#!/bin/bash
set -euo pipefail

# =============================================================================
# Task 33: RNA-seq Isoform-level Analysis
#
# DAG (depth 7, nested diamond):
# L0: PE reads + genome + GTF
# L1: fastp (trim)
# L2: STAR genomeGenerate (index)
# L3: STAR 2-pass align (pass 1 → collect junctions → pass 2)
# L4: samtools sort + index
#     ├────────────────────────────────────────┐
# L5: StringTie (isoform assembly)      featureCounts (gene counts)
#     │                                        │
# L6: gffcompare (vs reference GTF)      gene count summary
#     │
# L7: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${REF}/genome.fa"
GTF="${REF}/genes.gtf"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,star_index,aligned,stringtie,gffcompare,counts} "${RES}"

# L1: Trim
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
fi

# L2: STAR index
log_step "L2: STAR genomeGenerate"
if [ ! -f "${OUT}/star_index/SA" ]; then
    STAR --runMode genomeGenerate --genomeDir "${OUT}/star_index" \
         --genomeFastaFiles "${GENOME}" --sjdbGTFfile "${GTF}" \
         --runThreadN ${THREADS} --genomeSAindexNbases 11
fi

# L3: STAR 2-pass align
log_step "L3: STAR align"
if [ ! -f "${OUT}/aligned/Aligned.sortedByCoord.out.bam" ]; then
    STAR --genomeDir "${OUT}/star_index" \
         --readFilesIn "${OUT}/trimmed/R1.fastq.gz" "${OUT}/trimmed/R2.fastq.gz" \
         --readFilesCommand zcat --runThreadN ${THREADS} \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${OUT}/aligned/" \
         --twopassMode Basic \
         --quantMode GeneCounts \
         --outSAMstrandField intronMotif
fi

# L4: Index BAM
log_step "L4: samtools index"
if [ ! -f "${OUT}/aligned/Aligned.sortedByCoord.out.bam.bai" ]; then
    samtools index "${OUT}/aligned/Aligned.sortedByCoord.out.bam"
fi
BAM="${OUT}/aligned/Aligned.sortedByCoord.out.bam"

# L5 LEFT: StringTie isoform assembly
log_step "L5-LEFT: StringTie"
if [ ! -f "${OUT}/stringtie/assembled.gtf" ]; then
    stringtie "${BAM}" -G "${GTF}" -o "${OUT}/stringtie/assembled.gtf" \
              -p ${THREADS} -e -A "${OUT}/stringtie/gene_abund.tab"
fi

# L5 RIGHT: featureCounts gene-level
log_step "L5-RIGHT: featureCounts"
if [ ! -f "${OUT}/counts/gene_counts.txt" ]; then
    featureCounts -a "${GTF}" -o "${OUT}/counts/gene_counts.txt" \
                  -p --countReadPairs -T ${THREADS} "${BAM}"
fi

# L6: gffcompare
log_step "L6: gffcompare"
if [ ! -f "${OUT}/gffcompare/gffcmp.stats" ]; then
    gffcompare -r "${GTF}" -o "${OUT}/gffcompare/gffcmp" "${OUT}/stringtie/assembled.gtf"
fi

# L7: MERGE
log_step "L7-MERGE"

# STAR stats
TOTAL_READS=$(grep "Number of input reads" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")
MAPPED_READS=$(grep "Uniquely mapped reads number" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")
MAPPED_PCT=$(grep "Uniquely mapped reads %" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' | tr -d '%' || echo "0")
MULTI_PCT=$(grep "% of reads mapped to multiple loci" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' | tr -d '%' || echo "0")
SPLICE_JUNCTIONS=$(grep "Number of splices: Total" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")

# StringTie stats
TRANSCRIPTS=$(grep -c "transcript" "${OUT}/stringtie/assembled.gtf" 2>/dev/null || true)
TRANSCRIPTS=${TRANSCRIPTS:-0}
GENES_EXPRESSED=$(awk -F'\t' 'NR>1 && $9 > 0' "${OUT}/stringtie/gene_abund.tab" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# gffcompare stats
SENSITIVITY=$(grep "Sensitivity" "${OUT}/gffcompare/gffcmp.stats" 2>/dev/null | head -1 | awk '{print $3}' || echo "0")
PRECISION=$(grep "Precision" "${OUT}/gffcompare/gffcmp.stats" 2>/dev/null | head -1 | awk '{print $5}' || echo "0")

# featureCounts stats
ASSIGNED=$(grep "Assigned" "${OUT}/counts/gene_counts.txt.summary" | awk '{print $2}' || echo "0")
GENES_WITH_COUNTS=$(awk 'NR>2 && $NF>0' "${OUT}/counts/gene_counts.txt" | wc -l | tr -d ' ')

cat > "${RES}/isoform_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
uniquely_mapped,${MAPPED_READS}
unique_mapping_pct,${MAPPED_PCT}
multimapped_pct,${MULTI_PCT}
splice_junctions,${SPLICE_JUNCTIONS}
assembled_transcripts,${TRANSCRIPTS}
genes_expressed,${GENES_EXPRESSED}
transcript_sensitivity,${SENSITIVITY}
transcript_precision,${PRECISION}
assigned_read_pairs,${ASSIGNED}
genes_with_counts,${GENES_WITH_COUNTS}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/isoform_report.csv"
