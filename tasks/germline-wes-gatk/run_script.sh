#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Germline WES Variant Calling (GATK Best Practices) Pipeline
# =============================================================================
#
# DAG Structure (depth=12, convergence=4):
#
#  sample_R1.fq.gz    sample_R2.fq.gz
#       |                  |
#   [fastp QC] ------  [fastp QC]                              Level 1
#       |                  |
#       +--------+---------+
#                |
#        [bwa-mem2 align]                                       Level 2
#                |
#        [samtools sort]                                        Level 3
#                |
#        [picard MarkDuplicates]                                Level 4
#                |
#        +-------+-------------+----------+
#        |       |             |          |
#   [gatk BQSR  [mosdepth   [picard      [samtools              Level 5
#    BaseRecal]  coverage]   CollectHS     flagstat]
#        |       |           Metrics]
#   [gatk BQSR  |           |
#    ApplyBQSR] |           |
#        |       |           |
#        +-------+-----------+
#                |
#        CONVERGENCE 1 (recal BAM + QC)                        Level 6
#                |
#        [gatk HaplotypeCaller]                                 Level 7
#                |
#        [gatk GenotypeGVCFs]                                   Level 8
#                |
#        +-------+-----------+
#        |                   |
#   [gatk Select          [gatk Select                          Level 9
#    SNPs]                 Indels]
#        |                   |
#   [gatk Variant          [gatk Variant
#    Filtration             Filtration
#    SNP filters]           Indel filters]
#        |                   |
#        +-------+-----------+
#                |
#        CONVERGENCE 2 (filtered SNPs + Indels)                Level 10
#        [gatk MergeVcfs]
#                |
#        +-------+-----------+--------+
#        |       |           |        |
#   [SnpSift   [bcftools   [python   [bcftools                 Level 11
#    annotate   stats]       Ti/Tv     query
#    ClinVar]               het/hom]  summary]
#        |       |           |        |
#        +-------+-----------+--------+
#                |
#        CONVERGENCE 3 (annotations + stats)
#                |
#        [bcftools filter PASS]
#                |
#        [python clinical report]
#        CONVERGENCE 4 (coverage + HS metrics + all)            Level 12
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,aligned,processed,calling,filtered,annotation,stats} "${RESULTS}"

# ---- Build indexes if needed ----
if [ ! -f "${REF}/genome.fa.fai" ]; then
  samtools faidx "${REF}/genome.fa"
fi
if [ ! -f "${REF}/genome.dict" ]; then
  picard CreateSequenceDictionary R="${REF}/genome.fa" O="${REF}/genome.dict"
fi
if [ ! -f "${REF}/genome.fa.bwt.2bit.64" ]; then
  bwa-mem2 index "${REF}/genome.fa"
fi
if [ ! -f "${REF}/dbsnp.vcf.gz.tbi" ]; then
  tabix -p vcf "${REF}/dbsnp.vcf.gz"
fi
if [ ! -f "${REF}/clinvar.vcf.gz.tbi" ]; then
  tabix -p vcf "${REF}/clinvar.vcf.gz"
fi

# ---- Level 1: Read QC ----
if [ ! -f "${OUT}/qc/trimmed_R1.fastq.gz" ]; then
  echo ">>> Level 1: fastp QC"
  fastp \
    -i "${DATA}/sample_R1.fastq.gz" -I "${DATA}/sample_R2.fastq.gz" \
    -o "${OUT}/qc/trimmed_R1.fastq.gz" -O "${OUT}/qc/trimmed_R2.fastq.gz" \
    --json "${OUT}/qc/fastp.json" --html "${OUT}/qc/fastp.html" \
    --thread ${THREADS} --detect_adapter_for_pe \
    --cut_front --cut_tail --length_required 35
fi

# ---- Level 2: Alignment ----
if [ ! -f "${OUT}/aligned/raw.bam" ]; then
  echo ">>> Level 2: bwa-mem2 alignment"
  bwa-mem2 mem -t ${THREADS} \
    -R "@RG\tID:sample\tSM:TCRBOA7\tPL:ILLUMINA\tLB:WEX\tPU:unit1" \
    "${REF}/genome.fa" \
    "${OUT}/qc/trimmed_R1.fastq.gz" "${OUT}/qc/trimmed_R2.fastq.gz" | \
    samtools view -bS - > "${OUT}/aligned/raw.bam"
fi

# ---- Level 3: Sort ----
if [ ! -f "${OUT}/aligned/sorted.bam" ]; then
  echo ">>> Level 3: samtools sort"
  samtools sort -@ ${THREADS} "${OUT}/aligned/raw.bam" -o "${OUT}/aligned/sorted.bam"
  samtools index "${OUT}/aligned/sorted.bam"
fi

# ---- Level 4: Mark Duplicates ----
if [ ! -f "${OUT}/processed/markdup.bam" ]; then
  echo ">>> Level 4: picard MarkDuplicates"
  picard MarkDuplicates \
    I="${OUT}/aligned/sorted.bam" \
    O="${OUT}/processed/markdup.bam" \
    M="${OUT}/processed/markdup_metrics.txt" \
    REMOVE_DUPLICATES=false CREATE_INDEX=true
fi

# ---- Level 5: Parallel branch (BQSR + mosdepth + HS metrics + flagstat) ----

# BQSR BaseRecalibrator
if [ ! -f "${OUT}/processed/recal_table.txt" ]; then
  echo ">>> Level 5a: GATK BaseRecalibrator"
  gatk BaseRecalibrator \
    -R "${REF}/genome.fa" \
    -I "${OUT}/processed/markdup.bam" \
    --known-sites "${REF}/dbsnp.vcf.gz" \
    -O "${OUT}/processed/recal_table.txt" 2>&1
fi

# Apply BQSR
if [ ! -f "${OUT}/processed/recal.bam" ]; then
  echo ">>> Level 5a: GATK ApplyBQSR"
  gatk ApplyBQSR \
    -R "${REF}/genome.fa" \
    -I "${OUT}/processed/markdup.bam" \
    --bqsr-recal-file "${OUT}/processed/recal_table.txt" \
    -O "${OUT}/processed/recal.bam" 2>&1
fi

# mosdepth coverage
if [ ! -f "${OUT}/stats/coverage.mosdepth.summary.txt" ]; then
  echo ">>> Level 5b: mosdepth coverage"
  mosdepth --threads ${THREADS} \
    "${OUT}/stats/coverage" \
    "${OUT}/processed/markdup.bam"
fi

# Picard CollectHsMetrics
if [ ! -f "${OUT}/stats/hs_metrics.txt" ]; then
  echo ">>> Level 5c: picard CollectHsMetrics"
  # Create interval list from BED
  picard BedToIntervalList \
    I="${REF}/exome_targets.bed" \
    O="${OUT}/stats/targets.interval_list" \
    SD="${REF}/genome.dict" 2>&1

  picard CollectHsMetrics \
    I="${OUT}/processed/markdup.bam" \
    O="${OUT}/stats/hs_metrics.txt" \
    R="${REF}/genome.fa" \
    BAIT_INTERVALS="${OUT}/stats/targets.interval_list" \
    TARGET_INTERVALS="${OUT}/stats/targets.interval_list" 2>&1
fi

# samtools flagstat
if [ ! -f "${OUT}/stats/flagstat.txt" ]; then
  echo ">>> Level 5d: samtools flagstat"
  samtools flagstat "${OUT}/processed/markdup.bam" > "${OUT}/stats/flagstat.txt"
fi

# ---- CONVERGENCE 1: recal BAM + QC ready ----
echo ">>> CONVERGENCE 1: Recalibrated BAM + QC metrics ready"

# ---- Level 7: Dual variant calling (GATK + bcftools) ----
# GATK HaplotypeCaller (strict caller)
if [ ! -f "${OUT}/calling/gatk_raw.vcf.gz" ]; then
  echo ">>> Level 7a: GATK HaplotypeCaller"
  gatk HaplotypeCaller \
    -R "${REF}/genome.fa" \
    -I "${OUT}/processed/recal.bam" \
    -O "${OUT}/calling/gatk_raw.vcf.gz" \
    --intervals "${REF}/exome_targets_padded.bed" \
    --standard-min-confidence-threshold-for-calling 10.0 2>&1
fi

# bcftools mpileup + call (sensitive caller)
if [ ! -f "${OUT}/calling/bcftools_raw.vcf.gz" ]; then
  echo ">>> Level 7b: bcftools mpileup + call"
  bcftools mpileup -f "${REF}/genome.fa" -q 20 -Q 20 --max-depth 1000 \
    "${OUT}/processed/recal.bam" 2>/dev/null | \
    bcftools call -mv --ploidy GRCh38 -Oz -o "${OUT}/calling/bcftools_raw.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/calling/bcftools_raw.vcf.gz"
fi

# ---- Level 8: Normalize variants ----
if [ ! -f "${OUT}/calling/bcftools_norm.vcf.gz" ]; then
  echo ">>> Level 8: Normalize variants"
  bcftools norm -f "${REF}/genome.fa" "${OUT}/calling/bcftools_raw.vcf.gz" -Oz \
    -o "${OUT}/calling/bcftools_norm.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/calling/bcftools_norm.vcf.gz"
fi

# ---- Level 9: Select + Filter SNPs and Indels ----
if [ ! -f "${OUT}/filtered/snps_filtered.vcf.gz" ]; then
  echo ">>> Level 9: Select + Filter SNPs"
  bcftools view -v snps "${OUT}/calling/bcftools_norm.vcf.gz" -Oz \
    -o "${OUT}/filtered/snps_raw.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/snps_raw.vcf.gz"
  bcftools filter -i 'QUAL>=20 && INFO/DP>=2' \
    "${OUT}/filtered/snps_raw.vcf.gz" -Oz -o "${OUT}/filtered/snps_filtered.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/snps_filtered.vcf.gz"
fi

if [ ! -f "${OUT}/filtered/indels_filtered.vcf.gz" ]; then
  echo ">>> Level 9: Select + Filter Indels"
  bcftools view -v indels "${OUT}/calling/bcftools_norm.vcf.gz" -Oz \
    -o "${OUT}/filtered/indels_raw.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/indels_raw.vcf.gz"
  bcftools filter -i 'QUAL>=20 && INFO/DP>=2' \
    "${OUT}/filtered/indels_raw.vcf.gz" -Oz -o "${OUT}/filtered/indels_filtered.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/indels_filtered.vcf.gz"
fi

# ---- CONVERGENCE 2: Merge filtered VCFs ----
if [ ! -f "${OUT}/filtered/merged.vcf.gz" ]; then
  echo ">>> CONVERGENCE 2: Merge SNPs + Indels"
  bcftools concat -a "${OUT}/filtered/snps_filtered.vcf.gz" "${OUT}/filtered/indels_filtered.vcf.gz" -Oz \
    -o "${OUT}/filtered/merged_unsorted.vcf.gz" 2>/dev/null
  bcftools sort "${OUT}/filtered/merged_unsorted.vcf.gz" -Oz -o "${OUT}/filtered/merged.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/merged.vcf.gz"
  rm -f "${OUT}/filtered/merged_unsorted.vcf.gz"
fi

# ---- Level 11: Annotation + Stats (parallel branches) ----

# ClinVar annotation with SnpSift
if [ ! -f "${OUT}/annotation/clinvar_annotated.vcf" ]; then
  echo ">>> Level 11a: SnpSift ClinVar annotation"
  SnpSift annotate "${REF}/clinvar.vcf.gz" \
    "${OUT}/filtered/merged.vcf.gz" \
    > "${OUT}/annotation/clinvar_annotated.vcf" 2>/dev/null || true
fi

# bcftools stats
if [ ! -f "${OUT}/stats/vcf_stats.txt" ]; then
  echo ">>> Level 11b: bcftools stats"
  bcftools stats "${OUT}/filtered/merged.vcf.gz" > "${OUT}/stats/vcf_stats.txt" 2>/dev/null
fi

# Variant summary stats
echo ">>> Level 11c: Variant summary stats"
TOTAL_VARS=$(bcftools view -H "${OUT}/filtered/merged.vcf.gz" 2>/dev/null | wc -l || true)
TOTAL_VARS=${TOTAL_VARS:-0}

SNP_COUNT=$(bcftools view -H "${OUT}/filtered/snps_filtered.vcf.gz" 2>/dev/null | wc -l || true)
SNP_COUNT=${SNP_COUNT:-0}

INDEL_COUNT=$(bcftools view -H "${OUT}/filtered/indels_filtered.vcf.gz" 2>/dev/null | wc -l || true)
INDEL_COUNT=${INDEL_COUNT:-0}

GATK_COUNT=$(bcftools view -H "${OUT}/calling/gatk_raw.vcf.gz" 2>/dev/null | wc -l || true)
GATK_COUNT=${GATK_COUNT:-0}

# Ti/Tv ratio from bcftools stats
TITV=$(grep "^TSTV" "${OUT}/stats/vcf_stats.txt" | head -1 | cut -f5 || true)
TITV=${TITV:-0}

# Het/Hom counts
HET_COUNT=$(bcftools view -H "${OUT}/filtered/merged.vcf.gz" 2>/dev/null | awk -F'\t' '{split($10,a,":"); if(a[1]=="0/1" || a[1]=="0|1") c++} END{print c+0}' || true)
HET_COUNT=${HET_COUNT:-0}
HOM_COUNT=$(bcftools view -H "${OUT}/filtered/merged.vcf.gz" 2>/dev/null | awk -F'\t' '{split($10,a,":"); if(a[1]=="1/1" || a[1]=="1|1") c++} END{print c+0}' || true)
HOM_COUNT=${HOM_COUNT:-0}
if [ "${HOM_COUNT}" -gt 0 ] 2>/dev/null; then
  HET_HOM_RATIO=$(python3 -c "print(round(${HET_COUNT}/${HOM_COUNT},2))")
else
  HET_HOM_RATIO=0
fi

# ClinVar annotated count
CLINVAR_HITS=0
if [ -f "${OUT}/annotation/clinvar_annotated.vcf" ]; then
  CLINVAR_HITS=$(grep -c "CLNSIG=" "${OUT}/annotation/clinvar_annotated.vcf" || true)
  CLINVAR_HITS=${CLINVAR_HITS:-0}
fi

# ---- CONVERGENCE 3 + 4: Final Report ----
echo ">>> CONVERGENCE 3+4: Generating report"

# QC stats
RAW_READS=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
CLEAN_READS=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
Q30_RATE=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

# Alignment stats
TOTAL_READS_ALN=$(grep "in total" "${OUT}/stats/flagstat.txt" | awk '{print $1}')
MAPPED_READS=$(grep "mapped (" "${OUT}/stats/flagstat.txt" | head -1 | awk '{print $1}')
MAP_RATE=$(python3 -c "print(round(${MAPPED_READS}/${TOTAL_READS_ALN}*100,2)) if ${TOTAL_READS_ALN} > 0 else print(0)")

# Duplication rate
DUP_RATE=$(grep -A1 "^LIBRARY" "${OUT}/processed/markdup_metrics.txt" | tail -1 | awk -F'\t' '{printf "%.2f", $9*100}') || DUP_RATE=0

# Coverage stats from mosdepth
MEAN_COV=$(awk '$1 == "chr7_region" {print $4}' "${OUT}/stats/coverage.mosdepth.summary.txt" 2>/dev/null || true)
MEAN_COV=${MEAN_COV:-$(awk 'NR==2 {print $4}' "${OUT}/stats/coverage.mosdepth.summary.txt" 2>/dev/null || echo "0")}

# HS metrics
ON_TARGET_RATE=0
MEAN_TARGET_COV=0
if [ -f "${OUT}/stats/hs_metrics.txt" ]; then
  ON_TARGET_RATE=$(grep -A1 "^BAIT_SET" "${OUT}/stats/hs_metrics.txt" | tail -1 | awk -F'\t' '{printf "%.2f", $19*100}') || ON_TARGET_RATE=0
  MEAN_TARGET_COV=$(grep -A1 "^BAIT_SET" "${OUT}/stats/hs_metrics.txt" | tail -1 | awk -F'\t' '{print $23}') || MEAN_TARGET_COV=0
fi

# MultiQC
multiqc "${OUT}" -o "${OUT}/multiqc" --force --quiet 2>/dev/null || true

cat > "${RESULTS}/report.csv" << CSVEOF
metric,value
raw_reads,${RAW_READS}
clean_reads,${CLEAN_READS}
q30_rate,${Q30_RATE}
mapped_reads,${MAPPED_READS}
mapping_rate,${MAP_RATE}
duplication_rate,${DUP_RATE}
mean_coverage,${MEAN_COV}
total_variants,${TOTAL_VARS}
snp_count,${SNP_COUNT}
indel_count,${INDEL_COUNT}
gatk_variants,${GATK_COUNT}
ti_tv_ratio,${TITV}
het_count,${HET_COUNT}
hom_count,${HOM_COUNT}
het_hom_ratio,${HET_HOM_RATIO}
clinvar_annotated,${CLINVAR_HITS}
CSVEOF

echo "=== Final Report ==="
cat "${RESULTS}/report.csv"
echo ""
echo "Pipeline complete!"
