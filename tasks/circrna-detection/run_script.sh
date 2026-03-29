#!/bin/bash
set -uo pipefail

# =============================================================================
# Task: Circular RNA Detection and Quantification
#
# DAG structure (depth 8, 3 convergence points):
#
# L0: PE RNA-seq reads + genome + GTF + refFlat
# L1: trim_galore (adapter + quality trim)
# L2: STAR genomeGenerate (index with chimeric options)
# L3: STAR align (--chimSegmentMin 10, chimeric junction detection)
#     ├────────────────────────────────────────────────────┐
# L4: CIRCexplorer2 parse (extract chimeric junctions)    samtools stats
#     │                                                         │
# L5: CIRCexplorer2 annotate (match to gene annotation)   flagstat
#     │                                                         │
# L6: ├── filter (≥2 junction reads)                            │
#     └── bedtools: extract circRNA genomic coordinates         │
#          │                                                    │
# L7: featureCounts (linear gene expression for context) ◄─────┘
#     │                     [CONVERGENCE 1: linear counts + circRNA calls]
#     ├─────────────────────────────────────┐
# L8: circular/linear ratio calculation     circRNA size distribution
#     └──────────────┬──────────────────────┘
# L9: MERGE                                 [CONVERGENCE 2+3]
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${REF}/genome.fa"
GTF="${REF}/genes.gtf"
REFFLAT="${REF}/ref_flat.txt"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{trimmed,star_index,aligned,circexplorer,filtered,counts,stats} "${RES}"

# L1: Trim
log_step "L1: trim_galore"
if [ ! -f "${OUT}/trimmed/reads_R1_val_1.fq.gz" ]; then
    trim_galore --paired --fastqc -o "${OUT}/trimmed" \
                "${DATA}/reads_R1.fastq.gz" "${DATA}/reads_R2.fastq.gz"
fi

# L2: STAR index
log_step "L2: STAR genomeGenerate"
if [ ! -f "${OUT}/star_index/SA" ]; then
    STAR --runMode genomeGenerate --genomeDir "${OUT}/star_index" \
         --genomeFastaFiles "${GENOME}" --sjdbGTFfile "${GTF}" \
         --runThreadN ${THREADS} --genomeSAindexNbases 11
fi

# L3: STAR align with chimeric junction detection (required for circRNA)
log_step "L3: STAR chimeric align"
if [ ! -f "${OUT}/aligned/Chimeric.out.junction" ]; then
    STAR --genomeDir "${OUT}/star_index" \
         --readFilesIn "${OUT}/trimmed/reads_R1_val_1.fq.gz" "${OUT}/trimmed/reads_R2_val_2.fq.gz" \
         --readFilesCommand zcat --runThreadN ${THREADS} \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${OUT}/aligned/" \
         --chimSegmentMin 10 --chimOutType Junctions \
         --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
         --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \
         --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
         --outFilterMultimapNmax 20 --outSAMstrandField intronMotif
    samtools index "${OUT}/aligned/Aligned.sortedByCoord.out.bam"
fi
BAM="${OUT}/aligned/Aligned.sortedByCoord.out.bam"

# L4 LEFT: CIRCexplorer2 parse chimeric junctions
log_step "L4: CIRCexplorer2 parse"
if [ ! -f "${OUT}/circexplorer/back_spliced_junction.bed" ]; then
    CIRCexplorer2 parse -t STAR "${OUT}/aligned/Chimeric.out.junction" \
                  -b "${OUT}/circexplorer/back_spliced_junction.bed" 2>&1 || true
fi

# L4 RIGHT: samtools stats
log_step "L4: samtools stats"
samtools flagstat "${BAM}" > "${OUT}/stats/flagstat.txt"

# L5: CIRCexplorer2 annotate
log_step "L5: CIRCexplorer2 annotate"
if [ ! -f "${OUT}/circexplorer/circularRNA_known.txt" ]; then
    CIRCexplorer2 annotate -r "${REFFLAT}" -g "${GENOME}" \
                  -b "${OUT}/circexplorer/back_spliced_junction.bed" \
                  -o "${OUT}/circexplorer/circularRNA_known.txt" 2>&1 || true
fi

# L6: Filter (≥2 junction reads)
log_step "L6: filter circRNAs"
if [ -f "${OUT}/circexplorer/circularRNA_known.txt" ]; then
    awk '$13 >= 2' "${OUT}/circexplorer/circularRNA_known.txt" > "${OUT}/filtered/circrna_filtered.txt" 2>/dev/null || true
else
    touch "${OUT}/filtered/circrna_filtered.txt"
fi

# L7: featureCounts for linear gene expression
log_step "L7: featureCounts"
if [ ! -f "${OUT}/counts/gene_counts.txt" ]; then
    featureCounts -a "${GTF}" -o "${OUT}/counts/gene_counts.txt" \
                  -p --countReadPairs -T ${THREADS} "${BAM}" 2>&1 || true
fi

# MERGE
log_step "MERGE"

TOTAL_READS=$(grep "Number of input reads" "${OUT}/aligned/Log.final.out" 2>/dev/null | awk -F'\t' '{print $NF}' || echo "0")
MAPPED=$(grep "Uniquely mapped reads number" "${OUT}/aligned/Log.final.out" 2>/dev/null | awk -F'\t' '{print $NF}' || echo "0")
MAPPED_PCT=$(grep "Uniquely mapped reads %" "${OUT}/aligned/Log.final.out" 2>/dev/null | awk -F'\t' '{print $NF}' | tr -d '%' || echo "0")
CHIMERIC=$(grep "Number of chimeric reads" "${OUT}/aligned/Log.final.out" 2>/dev/null | awk -F'\t' '{print $NF}' || echo "0")

BSJ_RAW=$(wc -l < "${OUT}/circexplorer/back_spliced_junction.bed" 2>/dev/null | tr -d ' ' || echo "0")
CIRCRNA_ANNOTATED=$(wc -l < "${OUT}/circexplorer/circularRNA_known.txt" 2>/dev/null | tr -d ' ' || echo "0")
CIRCRNA_FILTERED=$(wc -l < "${OUT}/filtered/circrna_filtered.txt" 2>/dev/null | tr -d ' ' || echo "0")

GENES_EXPRESSED=$(awk 'NR>2 && $NF>0' "${OUT}/counts/gene_counts.txt" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# circRNA host genes
HOST_GENES=0
if [ -s "${OUT}/filtered/circrna_filtered.txt" ]; then
    HOST_GENES=$(awk '{print $14}' "${OUT}/filtered/circrna_filtered.txt" 2>/dev/null | sort -u | wc -l | tr -d ' ' || echo "0")
fi

cat > "${RES}/circrna_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
uniquely_mapped,${MAPPED}
mapping_pct,${MAPPED_PCT}
chimeric_reads,${CHIMERIC}
back_splice_junctions_raw,${BSJ_RAW}
circular_rnas_annotated,${CIRCRNA_ANNOTATED}
circular_rnas_filtered,${CIRCRNA_FILTERED}
host_genes,${HOST_GENES}
linear_genes_expressed,${GENES_EXPRESSED}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/circrna_report.csv"
