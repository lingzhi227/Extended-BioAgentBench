#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Nascent Transcription (GRO-seq) Analysis Pipeline
# =============================================================================
#
# DAG Structure (depth=9, convergence=4):
#
#  SRX882903_T1.fq.gz  SRX882903_T2.fq.gz  SRX882904_T1.fq.gz  SRX882904_T2.fq.gz
#       |                    |                    |                    |
#       +--------+-----------+                    +---------+----------+
#                |                                          |
#        [cat merge CD4+]                           [cat merge Jurkat]
#                |                                          |
#        [fastp QC] ----------------------------- [fastp QC]           Level 1
#                |                                          |
#        [bowtie2 rRNA depletion] ----------- [bowtie2 rRNA depletion] Level 2
#                |                                          |
#        [STAR align genome] ---------------- [STAR align genome]      Level 3
#                |                                          |
#        [samtools sort + dedup + filter] --- [samtools sort+dedup]     Level 4
#                |                                          |
#        +-------+-------------------+                      |
#        |       |                   |                      |
#    [bedtools [bedtools          [samtools                  |          Level 5
#     genomecov  genomecov         flagstat]                 |
#     (+ strand)] (- strand)]                               |
#        |       |                   |                      |
#        | [CONVERGENCE 1]          |                      |
#        | (bedGraphToBigWig        |                      |
#        |  + and - strand)         |                      |
#        |       |                   |                      |
#    +---+---+   |                   |                      |
#    |       |   |                   |                      |
#  [HOMER  [HOMER  [deeptools                               |          Level 6
#   findPeaks findPeaks computeMatrix                       |
#   -style   -style    (TSS profile)]                       |
#   tss]     groseq]                                        |
#    |       |           |                                  |
#    +---+---+           |                                  |
#        |               |                                  |
#    [CONVERGENCE 2]     |                                  |          Level 7
#    (TSS + enhancer)    |                                  |
#        |               |                                  |
#    +---+-------+       |                                  |
#    |   |       |       |                                  |
#  [HOMER [bedtools [deeptools                              |          Level 8
#   annotatePeaks intersect plotHeatmap]                    |
#   .pl]  (gene     |                                      |
#         regions)] |                                      |
#    |   |       |       |                                  |
#    +---+-------+-------+                                  |
#        |                                                  |
#    [CONVERGENCE 3] (annotated TSS + enhancers + heatmaps) |
#        |                                                  |
#    [CONVERGENCE 4] <-- QC + flagstat + both conditions    |          Level 9
#    [python report]  <-------------------------------------+
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,rrna,aligned,filtered,coverage,peaks,heatmaps,multiqc} "${RESULTS}"

# ---- Build indexes if needed ----
if [ ! -f "${REF}/rRNA_index.1.bt2" ]; then
  echo ">>> Building rRNA bowtie2 index..."
  bowtie2-build "${REF}/human_rRNA.fa" "${REF}/rRNA_index"
fi

if [ ! -f "${REF}/star_index/Genome" ]; then
  echo ">>> Building STAR genome index..."
  mkdir -p "${REF}/star_index"
  STAR --runMode genomeGenerate \
    --genomeDir "${REF}/star_index" \
    --genomeFastaFiles "${REF}/GRCh38_chr21.fa" \
    --sjdbGTFfile "${REF}/genes_chr21.gtf" \
    --genomeSAindexNbases 11 \
    --runThreadN ${THREADS}
fi

# ---- Merge replicates per condition ----
echo ">>> Merging replicates per condition..."
if [ ! -f "${DATA}/cd4_merged.fastq.gz" ]; then
  cat "${DATA}/SRX882903_T1.fastq.gz" "${DATA}/SRX882903_T2.fastq.gz" > "${DATA}/cd4_merged.fastq.gz"
fi
if [ ! -f "${DATA}/jurkat_merged.fastq.gz" ]; then
  cat "${DATA}/SRX882904_T1.fastq.gz" "${DATA}/SRX882904_T2.fastq.gz" > "${DATA}/jurkat_merged.fastq.gz"
fi

# ==== Process each condition ====
for SAMPLE in cd4 jurkat; do
  echo "===== Processing ${SAMPLE} ====="
  INPUT="${DATA}/${SAMPLE}_merged.fastq.gz"

  # ---- Level 1: Read QC with fastp ----
  if [ ! -f "${OUT}/qc/${SAMPLE}_trimmed.fastq.gz" ]; then
    echo ">>> [${SAMPLE}] Level 1: fastp QC"
    fastp \
      -i "${INPUT}" \
      -o "${OUT}/qc/${SAMPLE}_trimmed.fastq.gz" \
      --json "${OUT}/qc/${SAMPLE}_fastp.json" \
      --html "${OUT}/qc/${SAMPLE}_fastp.html" \
      --thread ${THREADS} \
      --length_required 20 \
      --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20
  fi

  # ---- Level 2: rRNA depletion via bowtie2 ----
  if [ ! -f "${OUT}/rrna/${SAMPLE}_no_rrna.fastq.gz" ]; then
    echo ">>> [${SAMPLE}] Level 2: bowtie2 rRNA removal"
    bowtie2 \
      -x "${REF}/rRNA_index" \
      -U "${OUT}/qc/${SAMPLE}_trimmed.fastq.gz" \
      --un-gz "${OUT}/rrna/${SAMPLE}_no_rrna.fastq.gz" \
      -S /dev/null \
      --threads ${THREADS} \
      --very-sensitive 2> "${OUT}/rrna/${SAMPLE}_rrna_bowtie2.log"
    echo ">>> [${SAMPLE}] rRNA mapping stats:"
    tail -1 "${OUT}/rrna/${SAMPLE}_rrna_bowtie2.log"
  fi

  # ---- Level 3: Genome alignment with STAR ----
  if [ ! -f "${OUT}/aligned/${SAMPLE}_Aligned.sortedByCoord.out.bam" ]; then
    echo ">>> [${SAMPLE}] Level 3: STAR genome alignment"
    rm -rf "${OUT}/aligned/${SAMPLE}_STARtmp"
    STAR \
      --runThreadN ${THREADS} \
      --genomeDir "${REF}/star_index" \
      --readFilesIn "${OUT}/rrna/${SAMPLE}_no_rrna.fastq.gz" \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "${OUT}/aligned/${SAMPLE}_" \
      --outSAMattributes NH HI NM MD \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNmax 3 \
      --alignIntronMax 1 \
      --outTmpDir "${OUT}/aligned/${SAMPLE}_STARtmp"
    samtools index "${OUT}/aligned/${SAMPLE}_Aligned.sortedByCoord.out.bam"
  fi

  # ---- Level 4: Deduplication + filtering ----
  if [ ! -f "${OUT}/filtered/${SAMPLE}_dedup.bam" ]; then
    echo ">>> [${SAMPLE}] Level 4: Dedup + filter"
    # Mark duplicates (position-based for GRO-seq)
    samtools markdup -r \
      "${OUT}/aligned/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
      "${OUT}/filtered/${SAMPLE}_dedup_unsorted.bam"
    # Sort and index
    samtools sort -@ ${THREADS} \
      "${OUT}/filtered/${SAMPLE}_dedup_unsorted.bam" \
      -o "${OUT}/filtered/${SAMPLE}_dedup.bam"
    samtools index "${OUT}/filtered/${SAMPLE}_dedup.bam"
    rm -f "${OUT}/filtered/${SAMPLE}_dedup_unsorted.bam"
  fi

  # ---- Level 5: Flagstat ----
  if [ ! -f "${OUT}/filtered/${SAMPLE}_flagstat.txt" ]; then
    echo ">>> [${SAMPLE}] Level 5: samtools flagstat"
    samtools flagstat "${OUT}/filtered/${SAMPLE}_dedup.bam" > "${OUT}/filtered/${SAMPLE}_flagstat.txt"
  fi

done

# ==== Main analysis on CD4+ (primary condition) ====
MAIN_BAM="${OUT}/filtered/cd4_dedup.bam"

# ---- Level 5: Strand-specific coverage ----
echo ">>> Level 5: Strand-specific coverage"
if [ ! -f "${OUT}/coverage/cd4_plus.bedGraph" ]; then
  # Plus strand: for GRO-seq SE, -strand + gives forward strand coverage
  bedtools genomecov -ibam "${MAIN_BAM}" -bg -strand + \
    > "${OUT}/coverage/cd4_plus.bedGraph"
  # Minus strand
  bedtools genomecov -ibam "${MAIN_BAM}" -bg -strand - \
    > "${OUT}/coverage/cd4_minus.bedGraph"
  echo "Plus strand regions: $(wc -l < "${OUT}/coverage/cd4_plus.bedGraph")"
  echo "Minus strand regions: $(wc -l < "${OUT}/coverage/cd4_minus.bedGraph")"
fi

# ---- CONVERGENCE 1: bedGraphToBigWig ----
echo ">>> CONVERGENCE 1: bedGraphToBigWig"
for STRAND in plus minus; do
  if [ ! -f "${OUT}/coverage/cd4_${STRAND}.bw" ]; then
    # Sort bedGraph by chr, start
    sort -k1,1 -k2,2n "${OUT}/coverage/cd4_${STRAND}.bedGraph" \
      > "${OUT}/coverage/cd4_${STRAND}_sorted.bedGraph"
    bedGraphToBigWig \
      "${OUT}/coverage/cd4_${STRAND}_sorted.bedGraph" \
      "${REF}/chrom.sizes" \
      "${OUT}/coverage/cd4_${STRAND}.bw"
  fi
done

# ---- Level 6: HOMER peak calling ----
echo ">>> Level 6: HOMER setup"
# Create HOMER tag directory
if [ ! -d "${OUT}/peaks/cd4_tags" ]; then
  makeTagDirectory "${OUT}/peaks/cd4_tags" "${MAIN_BAM}" -format sam
fi

# TSS-style peaks (promoters/TSSs)
if [ ! -f "${OUT}/peaks/cd4_tss_peaks.txt" ]; then
  echo ">>> Level 6a: HOMER findPeaks -style tss"
  findPeaks "${OUT}/peaks/cd4_tags" -style tss \
    -o "${OUT}/peaks/cd4_tss_peaks.txt" 2>&1 || true
fi

# GRO-seq style peaks (transcription units/enhancers)
if [ ! -f "${OUT}/peaks/cd4_groseq_peaks.txt" ]; then
  echo ">>> Level 6b: HOMER findPeaks -style groseq"
  findPeaks "${OUT}/peaks/cd4_tags" -style groseq \
    -o "${OUT}/peaks/cd4_groseq_peaks.txt" 2>&1 || true
fi

# ---- Level 6c: deeptools TSS profile ----
if [ ! -f "${OUT}/heatmaps/tss_matrix.gz" ]; then
  echo ">>> Level 6c: deeptools computeMatrix"
  # Create a BED file of gene TSSs from GTF
  awk '$3 == "exon" && $7 == "+" {print $1"\t"$4-1"\t"$4"\t"$10"\t0\t"$7}' \
    "${REF}/genes_chr21.gtf" | tr -d '";' | sort -u | head -500 \
    > "${OUT}/heatmaps/tss_sites_plus.bed" || true
  awk '$3 == "exon" && $7 == "-" {print $1"\t"$5-1"\t"$5"\t"$10"\t0\t"$7}' \
    "${REF}/genes_chr21.gtf" | tr -d '";' | sort -u | head -500 \
    > "${OUT}/heatmaps/tss_sites_minus.bed" || true
  cat "${OUT}/heatmaps/tss_sites_plus.bed" "${OUT}/heatmaps/tss_sites_minus.bed" \
    | sort -k1,1 -k2,2n > "${OUT}/heatmaps/tss_sites.bed"

  if [ -s "${OUT}/heatmaps/tss_sites.bed" ]; then
    computeMatrix reference-point \
      --referencePoint TSS \
      -S "${OUT}/coverage/cd4_plus.bw" "${OUT}/coverage/cd4_minus.bw" \
      -R "${OUT}/heatmaps/tss_sites.bed" \
      --beforeRegionStartLength 2000 \
      --afterRegionStartLength 2000 \
      -o "${OUT}/heatmaps/tss_matrix.gz" \
      --skipZeros \
      -p ${THREADS} 2>&1 || true
  fi
fi

# ---- CONVERGENCE 2: TSS + enhancer peaks ----
echo ">>> CONVERGENCE 2: merging TSS + enhancer results"
TSS_COUNT=0
GROSEQ_COUNT=0
if [ -f "${OUT}/peaks/cd4_tss_peaks.txt" ]; then
  TSS_COUNT=$(grep -c "^chr" "${OUT}/peaks/cd4_tss_peaks.txt" || true)
  TSS_COUNT=${TSS_COUNT:-0}
fi
if [ -f "${OUT}/peaks/cd4_groseq_peaks.txt" ]; then
  GROSEQ_COUNT=$(grep -c "^chr" "${OUT}/peaks/cd4_groseq_peaks.txt" || true)
  GROSEQ_COUNT=${GROSEQ_COUNT:-0}
fi
echo "TSS peaks: ${TSS_COUNT}, GRO-seq peaks: ${GROSEQ_COUNT}"

# ---- Level 7: HOMER annotatePeaks ----
if [ ! -f "${OUT}/peaks/cd4_tss_annotated.txt" ] && [ -f "${OUT}/peaks/cd4_tss_peaks.txt" ]; then
  echo ">>> Level 7: HOMER annotatePeaks"
  annotatePeaks.pl "${OUT}/peaks/cd4_tss_peaks.txt" \
    "${REF}/GRCh38_chr21.fa" \
    -gtf "${REF}/genes_chr21.gtf" \
    > "${OUT}/peaks/cd4_tss_annotated.txt" 2>/dev/null || true
fi

# ---- Level 8a: Gene region intersection ----
echo ">>> Level 8: Gene region intersection"
# Convert GTF to gene BED
if [ ! -f "${OUT}/peaks/gene_regions.bed" ]; then
  awk '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$10"\t0\t"$7}' \
    "${REF}/genes_chr21.gtf" | tr -d '";' | sort -k1,1 -k2,2n | \
    bedtools merge -i - > "${OUT}/peaks/gene_regions.bed" || true
fi

# Intersect GRO-seq peaks with gene regions
if [ -f "${OUT}/peaks/cd4_groseq_peaks.txt" ] && [ -s "${OUT}/peaks/gene_regions.bed" ]; then
  # Convert HOMER peaks to BED
  grep "^chr" "${OUT}/peaks/cd4_groseq_peaks.txt" | \
    awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$5}' | \
    sort -k1,1 -k2,2n > "${OUT}/peaks/cd4_groseq.bed" || true

  if [ -s "${OUT}/peaks/cd4_groseq.bed" ]; then
    GENIC_PEAKS=$(bedtools intersect -a "${OUT}/peaks/cd4_groseq.bed" \
      -b "${OUT}/peaks/gene_regions.bed" -u | wc -l || true)
    INTERGENIC_PEAKS=$(bedtools intersect -a "${OUT}/peaks/cd4_groseq.bed" \
      -b "${OUT}/peaks/gene_regions.bed" -v | wc -l || true)
  else
    GENIC_PEAKS=0
    INTERGENIC_PEAKS=0
  fi
else
  GENIC_PEAKS=0
  INTERGENIC_PEAKS=0
fi
echo "Genic peaks: ${GENIC_PEAKS}, Intergenic peaks: ${INTERGENIC_PEAKS}"

# ---- Level 8b: deeptools plotHeatmap ----
if [ -f "${OUT}/heatmaps/tss_matrix.gz" ] && [ ! -f "${OUT}/heatmaps/tss_heatmap.png" ]; then
  echo ">>> Level 8b: deeptools plotHeatmap"
  plotHeatmap -m "${OUT}/heatmaps/tss_matrix.gz" \
    -o "${OUT}/heatmaps/tss_heatmap.png" \
    --colorMap RdBu_r \
    --whatToShow "heatmap and colorbar" 2>&1 || true
fi

# ---- CONVERGENCE 3 + 4: Final Report ----
echo ">>> CONVERGENCE 3+4: Generating final report"

# Gather alignment stats from STAR logs (flagstat only has mapped reads)
CD4_STAR_INPUT=$(grep "Number of input reads" "${OUT}/aligned/cd4_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t')
CD4_UNIQUE=$(grep "Uniquely mapped reads number" "${OUT}/aligned/cd4_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t')
CD4_MAP_PCT=$(grep "Uniquely mapped reads %" "${OUT}/aligned/cd4_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t%')
JURKAT_STAR_INPUT=$(grep "Number of input reads" "${OUT}/aligned/jurkat_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t')
JURKAT_UNIQUE=$(grep "Uniquely mapped reads number" "${OUT}/aligned/jurkat_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t')
JURKAT_MAP_PCT=$(grep "Uniquely mapped reads %" "${OUT}/aligned/jurkat_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t%')

# rRNA rates (strip % sign)
CD4_RRNA_RATE=$(grep "overall alignment rate" "${OUT}/rrna/cd4_rrna_bowtie2.log" | awk '{print $1}' | tr -d '%')
JURKAT_RRNA_RATE=$(grep "overall alignment rate" "${OUT}/rrna/jurkat_rrna_bowtie2.log" | awk '{print $1}' | tr -d '%')

# fastp stats
CD4_RAW=$(python3 -c "import json; d=json.load(open('${OUT}/qc/cd4_fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
CD4_CLEAN=$(python3 -c "import json; d=json.load(open('${OUT}/qc/cd4_fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
CD4_Q30=$(python3 -c "import json; d=json.load(open('${OUT}/qc/cd4_fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

JURKAT_RAW=$(python3 -c "import json; d=json.load(open('${OUT}/qc/jurkat_fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
JURKAT_CLEAN=$(python3 -c "import json; d=json.load(open('${OUT}/qc/jurkat_fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
JURKAT_Q30=$(python3 -c "import json; d=json.load(open('${OUT}/qc/jurkat_fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

# Coverage stats from BigWig
PLUS_REGIONS=$(wc -l < "${OUT}/coverage/cd4_plus.bedGraph")
MINUS_REGIONS=$(wc -l < "${OUT}/coverage/cd4_minus.bedGraph")

# Peak annotation stats
ANNOTATED_PEAKS=0
if [ -f "${OUT}/peaks/cd4_tss_annotated.txt" ]; then
  ANNOTATED_PEAKS=$(tail -n +2 "${OUT}/peaks/cd4_tss_annotated.txt" | wc -l || true)
  ANNOTATED_PEAKS=${ANNOTATED_PEAKS:-0}
fi

# Dedup read counts from flagstat
CD4_DEDUP=$(grep "in total" "${OUT}/filtered/cd4_flagstat.txt" | awk '{print $1}')
JURKAT_DEDUP=$(grep "in total" "${OUT}/filtered/jurkat_flagstat.txt" | awk '{print $1}')

# MultiQC
echo ">>> Running MultiQC"
multiqc "${OUT}" -o "${OUT}/multiqc" --force --quiet 2>/dev/null || true

# Write CSV report
cat > "${RESULTS}/report.csv" << CSVEOF
metric,value
cd4_raw_reads,${CD4_RAW}
cd4_clean_reads,${CD4_CLEAN}
cd4_q30_rate,${CD4_Q30}
cd4_rrna_rate,${CD4_RRNA_RATE}
cd4_non_rrna_reads,${CD4_STAR_INPUT}
cd4_uniquely_mapped_reads,${CD4_UNIQUE}
cd4_unique_mapping_rate,${CD4_MAP_PCT}
cd4_dedup_reads,${CD4_DEDUP}
jurkat_raw_reads,${JURKAT_RAW}
jurkat_clean_reads,${JURKAT_CLEAN}
jurkat_q30_rate,${JURKAT_Q30}
jurkat_rrna_rate,${JURKAT_RRNA_RATE}
jurkat_non_rrna_reads,${JURKAT_STAR_INPUT}
jurkat_uniquely_mapped_reads,${JURKAT_UNIQUE}
jurkat_unique_mapping_rate,${JURKAT_MAP_PCT}
jurkat_dedup_reads,${JURKAT_DEDUP}
tss_peaks,${TSS_COUNT}
transcription_unit_peaks,${GROSEQ_COUNT}
genic_peaks,${GENIC_PEAKS}
intergenic_peaks,${INTERGENIC_PEAKS}
plus_strand_coverage_regions,${PLUS_REGIONS}
minus_strand_coverage_regions,${MINUS_REGIONS}
annotated_peaks,${ANNOTATED_PEAKS}
CSVEOF

echo "=== Final Report ==="
cat "${RESULTS}/report.csv"
echo ""
echo "Pipeline complete!"
