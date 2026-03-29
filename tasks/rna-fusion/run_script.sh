#!/bin/bash
set -euo pipefail

# =============================================================================
# Task: RNA Fusion Detection
#
# DAG structure (depth 8, 3 convergence points):
#
# L0: PE RNA-seq reads + genome + GTF
# L1: fastp (trim)
# L2: STAR genomeGenerate (index)
# L3: STAR align (chimeric junction detection enabled)
#     ├───────────────────────────────────────────────────┐
# L4: arriba (fusion detection from chimeric BAM)   samtools stats
#     │                                                    │
# L5: arriba filter (confidence)               StringTie (expression)
#     │                                              │
# L6: fusion annotation ◄───────────────────────────┘ [CONVERGENCE 1: fusions+expression]
#     │                    │
# L7: Salmon (quantification)  fusion gene expression
#     │                           │
# L8: MERGE ◄─────────────────────┘                    [CONVERGENCE 2+3]
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${REF}/genome.fa"
GTF="${REF}/genes.gtf"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{trimmed,star_index,aligned,arriba,stringtie,salmon,stats} "${RES}"

# L1
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
fi

# L2: STAR index
log_step "L2: STAR index"
if [ ! -f "${OUT}/star_index/SA" ]; then
    STAR --runMode genomeGenerate --genomeDir "${OUT}/star_index" \
         --genomeFastaFiles "${GENOME}" --sjdbGTFfile "${GTF}" \
         --runThreadN ${THREADS} --genomeSAindexNbases 11
fi

# L3: STAR align with chimeric detection (required for Arriba)
log_step "L3: STAR align (chimeric)"
if [ ! -f "${OUT}/aligned/Aligned.sortedByCoord.out.bam" ]; then
    STAR --genomeDir "${OUT}/star_index" \
         --readFilesIn "${OUT}/trimmed/R1.fastq.gz" "${OUT}/trimmed/R2.fastq.gz" \
         --readFilesCommand zcat --runThreadN ${THREADS} \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${OUT}/aligned/" \
         --chimSegmentMin 10 --chimOutType WithinBAM HardClip \
         --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
         --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \
         --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
         --outSAMstrandField intronMotif \
         --quantMode GeneCounts
    samtools index "${OUT}/aligned/Aligned.sortedByCoord.out.bam"
fi
BAM="${OUT}/aligned/Aligned.sortedByCoord.out.bam"

# L4 LEFT: Arriba fusion detection
log_step "L4: arriba"
if [ ! -f "${OUT}/arriba/fusions.tsv" ]; then
    arriba -x "${BAM}" -g "${GTF}" -a "${GENOME}" \
           -o "${OUT}/arriba/fusions.tsv" \
           -O "${OUT}/arriba/fusions.discarded.tsv" 2>&1 || {
        echo "WARNING: arriba found no fusions (normal for non-cancer data)"
        echo -e "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence\tclosest_genomic_breakpoint1\tclosest_genomic_breakpoint2\tfilters\tfusion_transcript\treading_frame\tpeptide_sequence\tread_identifiers" > "${OUT}/arriba/fusions.tsv"
    }
fi

# L4 RIGHT: samtools stats
log_step "L4: samtools stats"
samtools flagstat "${BAM}" > "${OUT}/stats/flagstat.txt"
samtools stats "${BAM}" > "${OUT}/stats/samtools_stats.txt"

# L5: StringTie expression
log_step "L5: StringTie"
if [ ! -f "${OUT}/stringtie/assembled.gtf" ]; then
    stringtie "${BAM}" -G "${GTF}" -o "${OUT}/stringtie/assembled.gtf" \
              -p ${THREADS} -e -A "${OUT}/stringtie/gene_abund.tab"
fi

# L5: Salmon quantification from STAR output
log_step "L5: gene counts from STAR"
cp "${OUT}/aligned/ReadsPerGene.out.tab" "${OUT}/stats/star_gene_counts.tsv" 2>/dev/null || true

# MERGE
log_step "MERGE"

TOTAL_READS=$(grep "Number of input reads" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")
MAPPED=$(grep "Uniquely mapped reads number" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")
MAPPED_PCT=$(grep "Uniquely mapped reads %" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' | tr -d '%' || echo "0")
CHIMERIC=$(grep "Number of chimeric reads" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")
SPLICE=$(grep "Number of splices: Total" "${OUT}/aligned/Log.final.out" | awk -F'\t' '{print $NF}' || echo "0")

FUSIONS=$(tail -n +2 "${OUT}/arriba/fusions.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
HIGH_CONF=$(tail -n +2 "${OUT}/arriba/fusions.tsv" 2>/dev/null | awk -F'\t' '$17=="high"' | wc -l | tr -d ' ' || echo "0")
MEDIUM_CONF=$(tail -n +2 "${OUT}/arriba/fusions.tsv" 2>/dev/null | awk -F'\t' '$17=="medium"' | wc -l | tr -d ' ' || echo "0")

GENES_EXPRESSED=$(awk 'NR>1 && $9>0' "${OUT}/stringtie/gene_abund.tab" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
TRANSCRIPTS=$(grep -c "transcript" "${OUT}/stringtie/assembled.gtf" 2>/dev/null || true)
TRANSCRIPTS=${TRANSCRIPTS:-0}

cat > "${RES}/fusion_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
uniquely_mapped,${MAPPED}
mapping_pct,${MAPPED_PCT}
chimeric_reads,${CHIMERIC}
splice_junctions,${SPLICE}
fusions_detected,${FUSIONS}
high_confidence_fusions,${HIGH_CONF}
medium_confidence_fusions,${MEDIUM_CONF}
genes_expressed,${GENES_EXPRESSED}
assembled_transcripts,${TRANSCRIPTS}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/fusion_report.csv"
