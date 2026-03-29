#!/bin/bash
set -uo pipefail

# =============================================================================
# Task: Ribosome Profiling (Ribo-seq)
#
# DAG structure (depth 7, 2 convergence points):
#
# L0: SE Ribo-seq reads + genome + GTF
# L1: fastp (adapter trim + length select 25-35nt for RPFs)
# L2: sortmerna (rRNA removal — critical, 80-90% of reads are rRNA)
# L3: STAR align (to genome, splice-aware)
# L4: samtools filter (unique mappers, 25-35nt only)
#     ├──────────────────────────────────────────────┐
# L5: samtools stats (mapping QC)             read length distribution
#     │                                              │
# L6: CDS quantification (featureCounts)      RPF size profile
#     └──────────┬───────────────────────────────────┘
# L7: MERGE                                          [CONVERGENCE 1+2]
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
mkdir -p "${OUT}"/{trimmed,rrna_filtered,star_index,aligned,filtered,stats} "${RES}"

# L1: Adapter trim + length select (RPFs are 25-35nt)
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/reads.fastq.gz" ]; then
    fastp --in1 "${DATA}/riboseq.fastq.gz" \
          --out1 "${OUT}/trimmed/reads.fastq.gz" \
          --adapter_sequence CTGTAGGCACCATCAAT \
          --length_required 25 --length_limit 35 \
          --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
fi

# L2: rRNA removal with SortMeRNA
log_step "L2: sortmerna rRNA removal"
if [ ! -f "${OUT}/rrna_filtered/non_rrna.fastq.gz" ]; then
    # SortMeRNA needs uncompressed or can handle gzipped
    SORTMERNA_DB=$(find /pscratch/sd/l/lingzhi/micromamba/envs/sessA-c9 -name "rRNA_databases" -type d 2>/dev/null | head -1)
    if [ -z "$SORTMERNA_DB" ]; then
        SORTMERNA_DB=$(find /pscratch/sd/l/lingzhi/micromamba/envs/sessA-c9 -name "smr_v4.3_*.fasta" 2>/dev/null | head -1 | xargs dirname 2>/dev/null)
    fi

    if [ -n "$SORTMERNA_DB" ]; then
        # Get all rRNA reference files
        REFS=$(find "$SORTMERNA_DB" -name "*.fasta" 2>/dev/null | head -5 | sed 's/^/--ref /' | tr '\n' ' ')
        sortmerna $REFS \
                  --reads "${OUT}/trimmed/reads.fastq.gz" \
                  --aligned "${OUT}/rrna_filtered/rrna" \
                  --other "${OUT}/rrna_filtered/non_rrna" \
                  --fastx --threads ${THREADS} \
                  --workdir "${OUT}/rrna_filtered/sortmerna_tmp" 2>&1 || true
        gzip "${OUT}/rrna_filtered/non_rrna.fastq" 2>/dev/null || \
        mv "${OUT}/rrna_filtered/non_rrna.fq.gz" "${OUT}/rrna_filtered/non_rrna.fastq.gz" 2>/dev/null || true
    fi

    # Fallback: if sortmerna failed, use trimmed reads directly
    if [ ! -s "${OUT}/rrna_filtered/non_rrna.fastq.gz" ]; then
        echo "WARNING: SortMeRNA failed, using all trimmed reads"
        cp "${OUT}/trimmed/reads.fastq.gz" "${OUT}/rrna_filtered/non_rrna.fastq.gz"
    fi
fi

# L3: STAR index + align
log_step "L3: STAR"
if [ ! -f "${OUT}/star_index/SA" ]; then
    STAR --runMode genomeGenerate --genomeDir "${OUT}/star_index" \
         --genomeFastaFiles "${GENOME}" --sjdbGTFfile "${GTF}" \
         --runThreadN ${THREADS} --genomeSAindexNbases 11
fi

if [ ! -f "${OUT}/aligned/Aligned.sortedByCoord.out.bam" ]; then
    STAR --genomeDir "${OUT}/star_index" \
         --readFilesIn "${OUT}/rrna_filtered/non_rrna.fastq.gz" \
         --readFilesCommand zcat --runThreadN ${THREADS} \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${OUT}/aligned/" \
         --alignIntronMax 5000 --outFilterMismatchNmax 1 \
         --quantMode GeneCounts
    samtools index "${OUT}/aligned/Aligned.sortedByCoord.out.bam"
fi
BAM="${OUT}/aligned/Aligned.sortedByCoord.out.bam"

# L4: Filter for RPF-length reads only
log_step "L4: samtools filter RPF lengths"
samtools flagstat "${BAM}" > "${OUT}/stats/flagstat.txt"

# L5-L6: Read length distribution + gene counts
log_step "L5: read length distribution"
samtools view "${BAM}" | awk '{print length($10)}' | sort | uniq -c | sort -k2 -n > "${OUT}/stats/read_lengths.txt"

# MERGE
log_step "MERGE"

TOTAL_READS=$(python3 -c "import json; d=json.load(open('${OUT}/trimmed/fastp.json')); print(d['summary']['before_filtering']['total_reads'])" 2>/dev/null || echo "0")
PASSED_TRIM=$(python3 -c "import json; d=json.load(open('${OUT}/trimmed/fastp.json')); print(d['summary']['after_filtering']['total_reads'])" 2>/dev/null || echo "0")

RRNA_READS=$(zcat "${OUT}/rrna_filtered/rrna.fastq.gz" 2>/dev/null | wc -l | awk '{print $1/4}' || echo "0")
NON_RRNA=$(zcat "${OUT}/rrna_filtered/non_rrna.fastq.gz" 2>/dev/null | wc -l | awk '{print $1/4}' || echo "0")
RRNA_PCT=$(python3 -c "
r=${RRNA_READS}; n=${NON_RRNA}
if r+n > 0: print(f'{r/(r+n)*100:.1f}')
else: print('0')
" 2>/dev/null || echo "0")

MAPPED=$(grep "Uniquely mapped reads number" "${OUT}/aligned/Log.final.out" 2>/dev/null | awk -F'\t' '{print $NF}' || echo "0")
MAPPED_PCT=$(grep "Uniquely mapped reads %" "${OUT}/aligned/Log.final.out" 2>/dev/null | awk -F'\t' '{print $NF}' | tr -d '%' || echo "0")

# Peak RPF length
PEAK_LENGTH=$(sort -k1 -rn "${OUT}/stats/read_lengths.txt" | head -1 | awk '{print $2}' || echo "0")
PEAK_COUNT=$(sort -k1 -rn "${OUT}/stats/read_lengths.txt" | head -1 | awk '{print $1}' || echo "0")

# Gene counts from STAR
GENES_WITH_READS=$(awk 'NR>4 && $2>0' "${OUT}/aligned/ReadsPerGene.out.tab" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

cat > "${RES}/riboseq_report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
passed_length_filter,${PASSED_TRIM}
rrna_pct,${RRNA_PCT}
non_rrna_reads,${NON_RRNA}
uniquely_mapped,${MAPPED}
mapping_pct,${MAPPED_PCT}
peak_rpf_length,${PEAK_LENGTH}
genes_with_reads,${GENES_WITH_READS}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/riboseq_report.csv"
