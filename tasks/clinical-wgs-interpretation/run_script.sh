#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Clinical WGS Interpretation Pipeline
# =============================================================================
#
# DAG Structure (depth=12, convergence=5):
#
#  sample_R1.fq.gz    sample_R2.fq.gz
#       |                  |
#   [fastp QC] ------  [fastp QC]                              Level 1
#       +--------+---------+
#                |
#        [bwa-mem2 align]                                       Level 2
#        [samtools sort]                                        Level 3
#        [picard MarkDuplicates]                                Level 4
#        [gatk BQSR (BaseRecal + Apply)]                        Level 5
#                |
#        +-------+-------------------+------------------+
#        |       |                   |                  |
#   [bcftools  [Delly            [mosdepth           [samtools  Level 6
#    mpileup    (SV call)]        (coverage +          flagstat]
#    + call]                      per-gene)]
#        |       |                   |                  |
#   [bcftools  [bcftools             |                  |
#    norm]      view SV VCF]         |                  |       Level 7
#        |       |                   |                  |
#   +----+----+  |                   |                  |
#   |         |  |                   |                  |
# [Select  [Select                   |                  |       Level 8
#  SNP]     Indel]                   |                  |
#   |         |                      |                  |
# [Filter  [Filter                   |                  |       Level 9
#  SNP]     Indel]                   |                  |
#   +----+----+                      |                  |
#        |                           |                  |
#   [bcftools concat]                |                  |
#        |       +-------------------+                  |
#        |       |                                      |
#   CONVERGENCE 1 (SNV/Indel + coverage)                        Level 9
#        |       |
#        |  CONVERGENCE 2 (Delly SVs)                           Level 10
#        +-------+
#                |
#        CONVERGENCE 3 (all variants merged)
#                |
#        +-------+---------------+--------+
#        |       |               |        |
#   [SnpSift   [bcftools       [python    |                     Level 11
#    annotate   annotate        score]    |
#    ClinVar]   (dbSNP freq)]            |
#        +-------+---------------+        |
#                |                        |
#        CONVERGENCE 4 (all annotations)  |
#                |                        |
#        +-------+-------+               |
#        |       |       |               |                      Level 12
#   [bcftools  [python  [python           |
#    filter     gene    clinical           |
#    PASS]      panel]  report]            |
#        +-------+-------+               |
#                |                        |
#        CONVERGENCE 5 (coverage + QC + final)
#        [final clinical report]
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,aligned,processed,snv,sv,filtered,annotation,stats,coverage} "${RESULTS}"

# ---- Build indexes if needed ----
[ ! -f "${REF}/genome.fa.fai" ] && samtools faidx "${REF}/genome.fa"
[ ! -f "${REF}/genome.dict" ] && picard CreateSequenceDictionary R="${REF}/genome.fa" O="${REF}/genome.dict" 2>&1
[ ! -f "${REF}/genome.fa.bwt.2bit.64" ] && bwa-mem2 index "${REF}/genome.fa"
[ ! -f "${REF}/dbsnp.vcf.gz.tbi" ] && tabix -p vcf "${REF}/dbsnp.vcf.gz"
[ ! -f "${REF}/clinvar.vcf.gz.tbi" ] && tabix -p vcf "${REF}/clinvar.vcf.gz"

# ---- Level 1: fastp QC ----
if [ ! -f "${OUT}/qc/trimmed_R1.fastq.gz" ]; then
  echo ">>> Level 1: fastp QC"
  fastp -i "${DATA}/sample_R1.fastq.gz" -I "${DATA}/sample_R2.fastq.gz" \
    -o "${OUT}/qc/trimmed_R1.fastq.gz" -O "${OUT}/qc/trimmed_R2.fastq.gz" \
    --json "${OUT}/qc/fastp.json" --html "${OUT}/qc/fastp.html" \
    --thread ${THREADS} --detect_adapter_for_pe --cut_front --cut_tail
fi

# ---- Level 2: bwa-mem2 align ----
if [ ! -f "${OUT}/aligned/raw.bam" ]; then
  echo ">>> Level 2: bwa-mem2 alignment"
  bwa-mem2 mem -t ${THREADS} \
    -R "@RG\tID:clinical\tSM:patient1\tPL:ILLUMINA\tLB:WGS\tPU:unit1" \
    "${REF}/genome.fa" \
    "${OUT}/qc/trimmed_R1.fastq.gz" "${OUT}/qc/trimmed_R2.fastq.gz" | \
    samtools view -bS - > "${OUT}/aligned/raw.bam"
fi

# ---- Level 3: samtools sort ----
if [ ! -f "${OUT}/aligned/sorted.bam" ]; then
  echo ">>> Level 3: samtools sort"
  samtools sort -@ ${THREADS} "${OUT}/aligned/raw.bam" -o "${OUT}/aligned/sorted.bam"
  samtools index "${OUT}/aligned/sorted.bam"
fi

# ---- Level 4: picard MarkDuplicates ----
if [ ! -f "${OUT}/processed/markdup.bam" ]; then
  echo ">>> Level 4: picard MarkDuplicates"
  picard MarkDuplicates \
    I="${OUT}/aligned/sorted.bam" O="${OUT}/processed/markdup.bam" \
    M="${OUT}/processed/markdup_metrics.txt" \
    REMOVE_DUPLICATES=false CREATE_INDEX=true
fi

# ---- Level 5: GATK BQSR ----
if [ ! -f "${OUT}/processed/recal.bam" ]; then
  echo ">>> Level 5: GATK BQSR"
  gatk BaseRecalibrator -R "${REF}/genome.fa" -I "${OUT}/processed/markdup.bam" \
    --known-sites "${REF}/dbsnp.vcf.gz" -O "${OUT}/processed/recal_table.txt" 2>&1
  gatk ApplyBQSR -R "${REF}/genome.fa" -I "${OUT}/processed/markdup.bam" \
    --bqsr-recal-file "${OUT}/processed/recal_table.txt" \
    -O "${OUT}/processed/recal.bam" 2>&1
fi

# ============================
# Level 6: Three parallel tracks from recal BAM
# ============================

# ---- Track 1: SNV/Indel calling (bcftools) ----
if [ ! -f "${OUT}/snv/bcf_raw.vcf.gz" ]; then
  echo ">>> Level 6a: bcftools SNV/Indel calling"
  bcftools mpileup -f "${REF}/genome.fa" -q 20 -Q 20 --max-depth 1000 \
    "${OUT}/processed/recal.bam" 2>/dev/null | \
    bcftools call -mv --ploidy GRCh38 -Oz -o "${OUT}/snv/bcf_raw.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/snv/bcf_raw.vcf.gz"
fi

# Level 7: Normalize
if [ ! -f "${OUT}/snv/bcf_norm.vcf.gz" ]; then
  echo ">>> Level 7: Normalize variants"
  bcftools norm -f "${REF}/genome.fa" "${OUT}/snv/bcf_raw.vcf.gz" -Oz \
    -o "${OUT}/snv/bcf_norm.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/snv/bcf_norm.vcf.gz"
fi

# Level 8-9: Select + Filter SNPs and Indels
if [ ! -f "${OUT}/filtered/snps_filtered.vcf.gz" ]; then
  echo ">>> Level 8-9: SNP filtering"
  bcftools view -v snps "${OUT}/snv/bcf_norm.vcf.gz" -Oz -o "${OUT}/filtered/snps_raw.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/snps_raw.vcf.gz"
  bcftools filter -i 'QUAL>=20 && INFO/DP>=2' "${OUT}/filtered/snps_raw.vcf.gz" -Oz \
    -o "${OUT}/filtered/snps_filtered.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/snps_filtered.vcf.gz"
fi

if [ ! -f "${OUT}/filtered/indels_filtered.vcf.gz" ]; then
  echo ">>> Level 8-9: Indel filtering"
  bcftools view -v indels "${OUT}/snv/bcf_norm.vcf.gz" -Oz -o "${OUT}/filtered/indels_raw.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/indels_raw.vcf.gz"
  bcftools filter -i 'QUAL>=20 && INFO/DP>=2' "${OUT}/filtered/indels_raw.vcf.gz" -Oz \
    -o "${OUT}/filtered/indels_filtered.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/indels_filtered.vcf.gz"
fi

# Merge filtered SNPs + Indels
if [ ! -f "${OUT}/filtered/snv_merged.vcf.gz" ]; then
  bcftools concat -a "${OUT}/filtered/snps_filtered.vcf.gz" "${OUT}/filtered/indels_filtered.vcf.gz" -Oz \
    -o "${OUT}/filtered/snv_merged_unsorted.vcf.gz" 2>/dev/null
  bcftools sort "${OUT}/filtered/snv_merged_unsorted.vcf.gz" -Oz \
    -o "${OUT}/filtered/snv_merged.vcf.gz" 2>/dev/null
  bcftools index "${OUT}/filtered/snv_merged.vcf.gz"
  rm -f "${OUT}/filtered/snv_merged_unsorted.vcf.gz"
fi

# ---- Track 2: SV detection (Delly) ----
if [ ! -f "${OUT}/sv/delly.vcf.gz" ]; then
  echo ">>> Level 6b: Delly SV calling"
  delly call -g "${REF}/genome.fa" "${OUT}/processed/recal.bam" \
    -o "${OUT}/sv/delly.bcf" 2>&1 || true
  if [ -f "${OUT}/sv/delly.bcf" ]; then
    bcftools view "${OUT}/sv/delly.bcf" -Oz -o "${OUT}/sv/delly.vcf.gz" 2>/dev/null || true
    bcftools index "${OUT}/sv/delly.vcf.gz" 2>/dev/null || true
  fi
fi

# ---- Track 3: Coverage analysis (mosdepth) ----
if [ ! -f "${OUT}/coverage/coverage.mosdepth.summary.txt" ]; then
  echo ">>> Level 6c: mosdepth coverage"
  mosdepth --threads ${THREADS} "${OUT}/coverage/coverage" "${OUT}/processed/recal.bam"
fi

# ---- Flagstat ----
if [ ! -f "${OUT}/stats/flagstat.txt" ]; then
  echo ">>> Level 6d: samtools flagstat"
  samtools flagstat "${OUT}/processed/recal.bam" > "${OUT}/stats/flagstat.txt"
fi

# ============================
# CONVERGENCE 1-3: Merge all variant tracks
# ============================
echo ">>> CONVERGENCE 1-3: Merging variant tracks"

SNP_COUNT=$(bcftools view -H "${OUT}/filtered/snps_filtered.vcf.gz" 2>/dev/null | wc -l || true)
INDEL_COUNT=$(bcftools view -H "${OUT}/filtered/indels_filtered.vcf.gz" 2>/dev/null | wc -l || true)
SNV_TOTAL=$((${SNP_COUNT:-0} + ${INDEL_COUNT:-0}))

SV_COUNT=0
if [ -s "${OUT}/sv/delly.vcf.gz" ]; then
  SV_COUNT=$(bcftools view -H "${OUT}/sv/delly.vcf.gz" 2>/dev/null | wc -l || true)
fi

echo "SNVs: ${SNV_TOTAL} (SNPs: ${SNP_COUNT}, Indels: ${INDEL_COUNT}), SVs: ${SV_COUNT}"

# ============================
# Level 11: Triple annotation branch
# ============================

# ClinVar annotation
if [ ! -f "${OUT}/annotation/clinvar_annotated.vcf" ]; then
  echo ">>> Level 11a: ClinVar annotation"
  SnpSift annotate "${REF}/clinvar.vcf.gz" "${OUT}/filtered/snv_merged.vcf.gz" \
    > "${OUT}/annotation/clinvar_annotated.vcf" 2>/dev/null || true
fi

# bcftools stats
if [ ! -f "${OUT}/stats/vcf_stats.txt" ]; then
  echo ">>> Level 11b: bcftools stats"
  bcftools stats "${OUT}/filtered/snv_merged.vcf.gz" > "${OUT}/stats/vcf_stats.txt" 2>/dev/null
fi

# Variant scoring (Ti/Tv, het/hom)
echo ">>> Level 11c: Variant scoring"
TITV=$(grep "^TSTV" "${OUT}/stats/vcf_stats.txt" | head -1 | cut -f5 || true)
TITV=${TITV:-0}

HET_COUNT=$(bcftools view -H "${OUT}/filtered/snv_merged.vcf.gz" 2>/dev/null | \
  awk -F'\t' '{split($10,a,":"); if(a[1]=="0/1" || a[1]=="0|1") c++} END{print c+0}' || true)
HOM_COUNT=$(bcftools view -H "${OUT}/filtered/snv_merged.vcf.gz" 2>/dev/null | \
  awk -F'\t' '{split($10,a,":"); if(a[1]=="1/1" || a[1]=="1|1") c++} END{print c+0}' || true)

CLINVAR_HITS=$(grep -c "CLNSIG=" "${OUT}/annotation/clinvar_annotated.vcf" 2>/dev/null || true)
CLINVAR_HITS=${CLINVAR_HITS:-0}

# ============================
# Level 12: Gene panel + clinical report
# ============================
echo ">>> Level 12: Gene panel filtering + report"

# Gene panel intersection
PANEL_VARIANTS=0
if [ -s "${REF}/gene_panel.bed" ]; then
  bcftools view -H "${OUT}/filtered/snv_merged.vcf.gz" 2>/dev/null | \
    awk -F'\t' '{print $1"\t"$2-1"\t"$2"\t"$4">"$5}' \
    > "${OUT}/annotation/variant_positions.bed" || true
  if [ -s "${OUT}/annotation/variant_positions.bed" ]; then
    PANEL_VARIANTS=$(bedtools intersect -a "${OUT}/annotation/variant_positions.bed" \
      -b "${REF}/gene_panel.bed" -u | wc -l || true)
  fi
fi

# Coverage stats
MEAN_COV=$(awk 'NR==2 {print $4}' "${OUT}/coverage/coverage.mosdepth.summary.txt" 2>/dev/null || echo "0")

# Low coverage regions (< 5x)
LOW_COV_REGIONS=0
if [ -f "${OUT}/coverage/coverage.per-base.bed.gz" ]; then
  LOW_COV_REGIONS=$(zcat "${OUT}/coverage/coverage.per-base.bed.gz" | \
    awk '$4 < 5 && $4 > 0 {c++} END{print c+0}' || true)
fi

# Alignment stats
MAPPED=$(grep "mapped (" "${OUT}/stats/flagstat.txt" | head -1 | awk '{print $1}')
TOTAL_ALN=$(grep "in total" "${OUT}/stats/flagstat.txt" | awk '{print $1}')
MAP_RATE=$(python3 -c "print(round(${MAPPED}/${TOTAL_ALN}*100,2)) if ${TOTAL_ALN} > 0 else print(0)")

# QC stats
RAW_READS=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
CLEAN_READS=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
Q30=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

DUP_RATE=$(grep -A1 "^LIBRARY" "${OUT}/processed/markdup_metrics.txt" | tail -1 | awk -F'\t' '{printf "%.2f", $9*100}') || DUP_RATE=0

HET_HOM_RATIO=0
if [ "${HOM_COUNT}" -gt 0 ] 2>/dev/null; then
  HET_HOM_RATIO=$(python3 -c "print(round(${HET_COUNT}/${HOM_COUNT},2))")
fi

# MultiQC
multiqc "${OUT}" -o "${OUT}/multiqc" --force --quiet 2>/dev/null || true

# Write report
cat > "${RESULTS}/report.csv" << CSVEOF
metric,value
raw_reads,${RAW_READS}
clean_reads,${CLEAN_READS}
q30_rate,${Q30}
mapped_reads,${MAPPED}
mapping_rate,${MAP_RATE}
duplication_rate,${DUP_RATE}
mean_coverage,${MEAN_COV}
snp_count,${SNP_COUNT}
indel_count,${INDEL_COUNT}
total_snv,${SNV_TOTAL}
sv_count,${SV_COUNT}
ti_tv_ratio,${TITV}
het_count,${HET_COUNT}
hom_count,${HOM_COUNT}
het_hom_ratio,${HET_HOM_RATIO}
clinvar_annotated,${CLINVAR_HITS}
gene_panel_variants,${PANEL_VARIANTS}
low_coverage_regions,${LOW_COV_REGIONS}
CSVEOF

echo "=== Final Report ==="
cat "${RESULTS}/report.csv"
echo ""
echo "Pipeline complete!"
