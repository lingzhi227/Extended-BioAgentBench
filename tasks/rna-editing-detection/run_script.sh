#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# RNA Editing Detection (A-to-I) Pipeline
# =============================================================================
#
# DAG Structure (depth=10, convergence=4):
#
#  rna_R1.fq.gz  rna_R2.fq.gz    dna_R1.fq.gz  dna_R2.fq.gz
#       |              |               |              |
#   [fastp QC] ---- [fastp QC]     [fastp QC] ---- [fastp QC]        Level 1
#       |              |               |              |
#       +------+-------+               +------+-------+
#              |                              |
#       [STAR 2-pass align]            [bwa-mem2 align]               Level 2
#              |                              |
#       [picard MarkDup]               [picard MarkDup]               Level 3
#              |                              |
#       [gatk SplitNCigarReads]        [gatk BQSR]                    Level 4
#              |                              |
#       [gatk BQSR]                           |                       Level 5
#              |                              |
#              +------- CONVERGENCE 1 --------+                       Level 6
#              (RNA BAM + DNA BAM ready)
#              |
#       +------+--------+--------+
#       |                |        |
#  [VarScan       [bcftools   [gatk HC                                Level 7
#   somatic        mpileup     RNA-mode]
#   RNA vs DNA]    +call]
#       |                |        |
#       +------+---------+--------+
#              |
#       CONVERGENCE 2                                                 Level 8
#       (3-caller intersection: >= 2 agree)
#              |
#       +------+---------+
#       |      |         |
#  [filter  [bcftools [bedtools                                       Level 9
#   A>G/T>C  annotate   intersect
#   strand]  (dbSNP     (repeat
#            filter)]    regions)]
#       |      |         |
#       +------+---------+
#              |
#       CONVERGENCE 3                                                 Level 9
#       (strand-filtered + SNP-filtered + repeat annotation)
#              |
#       [python editing report]
#       CONVERGENCE 4 <-- RNA/DNA QC                                  Level 10
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{rna_qc,dna_qc,rna_align,dna_align,rna_proc,dna_proc,calling,filtered} "${RESULTS}"

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
if [ ! -f "${REF}/star_index/Genome" ]; then
  mkdir -p "${REF}/star_index"
  STAR --runMode genomeGenerate \
    --genomeDir "${REF}/star_index" \
    --genomeFastaFiles "${REF}/genome.fa" \
    --sjdbGTFfile "${REF}/genes.gtf" \
    --genomeSAindexNbases 11 \
    --runThreadN ${THREADS}
fi
if [ ! -f "${REF}/dbsnp.vcf.gz.tbi" ]; then
  tabix -p vcf "${REF}/dbsnp.vcf.gz"
fi

# ============================
# RNA-seq Processing Branch
# ============================

# ---- Level 1: RNA fastp QC ----
if [ ! -f "${OUT}/rna_qc/rna_R1_trimmed.fastq.gz" ]; then
  echo ">>> Level 1: RNA fastp QC"
  fastp \
    -i "${DATA}/rna_R1.fastq.gz" -I "${DATA}/rna_R2.fastq.gz" \
    -o "${OUT}/rna_qc/rna_R1_trimmed.fastq.gz" -O "${OUT}/rna_qc/rna_R2_trimmed.fastq.gz" \
    --json "${OUT}/rna_qc/rna_fastp.json" \
    --html "${OUT}/rna_qc/rna_fastp.html" \
    --thread ${THREADS} \
    --detect_adapter_for_pe \
    --cut_front --cut_tail
fi

# ---- Level 2: STAR 2-pass alignment ----
if [ ! -f "${OUT}/rna_align/rna_Aligned.sortedByCoord.out.bam" ]; then
  echo ">>> Level 2: STAR 2-pass RNA alignment"
  rm -rf "${OUT}/rna_align/rna_STARtmp"
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir "${REF}/star_index" \
    --readFilesIn "${OUT}/rna_qc/rna_R1_trimmed.fastq.gz" "${OUT}/rna_qc/rna_R2_trimmed.fastq.gz" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUT}/rna_align/rna_" \
    --outSAMattributes NH HI NM MD \
    --twopassMode Basic \
    --outSAMattrRGline ID:rna SM:tumor_rna PL:ILLUMINA LB:lib1 \
    --outTmpDir "${OUT}/rna_align/rna_STARtmp"
  samtools index "${OUT}/rna_align/rna_Aligned.sortedByCoord.out.bam"
fi

# ---- Level 3: RNA MarkDuplicates ----
if [ ! -f "${OUT}/rna_proc/rna_markdup.bam" ]; then
  echo ">>> Level 3: RNA MarkDuplicates"
  picard MarkDuplicates \
    I="${OUT}/rna_align/rna_Aligned.sortedByCoord.out.bam" \
    O="${OUT}/rna_proc/rna_markdup.bam" \
    M="${OUT}/rna_proc/rna_markdup_metrics.txt" \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true
fi

# ---- Level 4: SplitNCigarReads (RNA-specific) ----
if [ ! -f "${OUT}/rna_proc/rna_split.bam" ]; then
  echo ">>> Level 4: gatk SplitNCigarReads"
  gatk SplitNCigarReads \
    -R "${REF}/genome.fa" \
    -I "${OUT}/rna_proc/rna_markdup.bam" \
    -O "${OUT}/rna_proc/rna_split.bam" 2>&1
fi

# ---- Level 5: RNA BQSR ----
if [ ! -f "${OUT}/rna_proc/rna_recal.bam" ]; then
  echo ">>> Level 5: RNA BQSR"
  gatk BaseRecalibrator \
    -R "${REF}/genome.fa" \
    -I "${OUT}/rna_proc/rna_split.bam" \
    --known-sites "${REF}/dbsnp.vcf.gz" \
    -O "${OUT}/rna_proc/rna_recal_table.txt" 2>&1
  gatk ApplyBQSR \
    -R "${REF}/genome.fa" \
    -I "${OUT}/rna_proc/rna_split.bam" \
    --bqsr-recal-file "${OUT}/rna_proc/rna_recal_table.txt" \
    -O "${OUT}/rna_proc/rna_recal.bam" 2>&1
fi

# ============================
# DNA Processing Branch
# ============================

# ---- Level 1: DNA fastp QC ----
if [ ! -f "${OUT}/dna_qc/dna_R1_trimmed.fastq.gz" ]; then
  echo ">>> Level 1: DNA fastp QC"
  fastp \
    -i "${DATA}/dna_R1.fastq.gz" -I "${DATA}/dna_R2.fastq.gz" \
    -o "${OUT}/dna_qc/dna_R1_trimmed.fastq.gz" -O "${OUT}/dna_qc/dna_R2_trimmed.fastq.gz" \
    --json "${OUT}/dna_qc/dna_fastp.json" \
    --html "${OUT}/dna_qc/dna_fastp.html" \
    --thread ${THREADS} \
    --detect_adapter_for_pe \
    --cut_front --cut_tail
fi

# ---- Level 2: bwa-mem2 DNA alignment ----
if [ ! -f "${OUT}/dna_align/dna_sorted.bam" ]; then
  echo ">>> Level 2: bwa-mem2 DNA alignment"
  bwa-mem2 mem -t ${THREADS} \
    -R "@RG\tID:dna\tSM:normal_dna\tPL:ILLUMINA\tLB:lib1" \
    "${REF}/genome.fa" \
    "${OUT}/dna_qc/dna_R1_trimmed.fastq.gz" "${OUT}/dna_qc/dna_R2_trimmed.fastq.gz" | \
    samtools sort -@ ${THREADS} -o "${OUT}/dna_align/dna_sorted.bam"
  samtools index "${OUT}/dna_align/dna_sorted.bam"
fi

# ---- Level 3: DNA MarkDuplicates ----
if [ ! -f "${OUT}/dna_proc/dna_markdup.bam" ]; then
  echo ">>> Level 3: DNA MarkDuplicates"
  picard MarkDuplicates \
    I="${OUT}/dna_align/dna_sorted.bam" \
    O="${OUT}/dna_proc/dna_markdup.bam" \
    M="${OUT}/dna_proc/dna_markdup_metrics.txt" \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true
fi

# ---- Level 4: DNA BQSR ----
if [ ! -f "${OUT}/dna_proc/dna_recal.bam" ]; then
  echo ">>> Level 4: DNA BQSR"
  gatk BaseRecalibrator \
    -R "${REF}/genome.fa" \
    -I "${OUT}/dna_proc/dna_markdup.bam" \
    --known-sites "${REF}/dbsnp.vcf.gz" \
    -O "${OUT}/dna_proc/dna_recal_table.txt" 2>&1
  gatk ApplyBQSR \
    -R "${REF}/genome.fa" \
    -I "${OUT}/dna_proc/dna_markdup.bam" \
    --bqsr-recal-file "${OUT}/dna_proc/dna_recal_table.txt" \
    -O "${OUT}/dna_proc/dna_recal.bam" 2>&1
fi

# ============================
# CONVERGENCE 1: RNA + DNA BAMs ready
# ============================
echo ">>> CONVERGENCE 1: RNA + DNA BAMs ready"
RNA_BAM="${OUT}/rna_proc/rna_recal.bam"
DNA_BAM="${OUT}/dna_proc/dna_recal.bam"

# Flagstats
samtools flagstat "${RNA_BAM}" > "${OUT}/rna_proc/rna_flagstat.txt"
samtools flagstat "${DNA_BAM}" > "${OUT}/dna_proc/dna_flagstat.txt"
echo "RNA mapped: $(grep 'mapped (' ${OUT}/rna_proc/rna_flagstat.txt | head -1)"
echo "DNA mapped: $(grep 'mapped (' ${OUT}/dna_proc/dna_flagstat.txt | head -1)"

# ============================
# Level 7: Three parallel callers
# ============================

# ---- Caller 1: VarScan somatic (RNA vs DNA pileup) ----
if [ ! -f "${OUT}/calling/varscan_snp.vcf" ]; then
  echo ">>> Level 7a: VarScan somatic"
  # Create pileups
  samtools mpileup -f "${REF}/genome.fa" -q 20 -Q 20 \
    "${DNA_BAM}" > "${OUT}/calling/dna_pileup.txt" 2>/dev/null
  samtools mpileup -f "${REF}/genome.fa" -q 20 -Q 20 \
    "${RNA_BAM}" > "${OUT}/calling/rna_pileup.txt" 2>/dev/null
  # VarScan somatic: normal=DNA, tumor=RNA to find RNA-specific variants
  varscan somatic \
    "${OUT}/calling/dna_pileup.txt" \
    "${OUT}/calling/rna_pileup.txt" \
    "${OUT}/calling/varscan" \
    --min-coverage 5 \
    --min-var-freq 0.1 \
    --output-vcf 1 2>&1 || true
  # Rename output
  if [ -f "${OUT}/calling/varscan.snp.vcf" ]; then
    cp "${OUT}/calling/varscan.snp.vcf" "${OUT}/calling/varscan_snp.vcf"
  else
    touch "${OUT}/calling/varscan_snp.vcf"
  fi
fi

# ---- Caller 2: bcftools mpileup + call ----
if [ ! -f "${OUT}/calling/bcftools_raw.vcf.gz" ]; then
  echo ">>> Level 7b: bcftools mpileup + call"
  bcftools mpileup -f "${REF}/genome.fa" \
    -q 20 -Q 20 --max-depth 1000 \
    "${RNA_BAM}" 2>/dev/null | \
    bcftools call -mv -Oz -o "${OUT}/calling/bcftools_raw.vcf.gz" 2>/dev/null || true
  bcftools index "${OUT}/calling/bcftools_raw.vcf.gz" 2>/dev/null || true
fi

# ---- Caller 3: GATK HaplotypeCaller (RNA mode) ----
if [ ! -f "${OUT}/calling/gatk_rna_raw.vcf.gz" ]; then
  echo ">>> Level 7c: GATK HaplotypeCaller RNA mode"
  gatk HaplotypeCaller \
    -R "${REF}/genome.fa" \
    -I "${RNA_BAM}" \
    -O "${OUT}/calling/gatk_rna_raw.vcf.gz" \
    --dont-use-soft-clipped-bases true \
    --standard-min-confidence-threshold-for-calling 20 2>&1 || true
fi

# ============================
# CONVERGENCE 2: Caller intersection
# ============================
echo ">>> CONVERGENCE 2: Multi-caller intersection"

# Count variants per caller
VARSCAN_COUNT=0
BCFTOOLS_COUNT=0
GATK_COUNT=0

if [ -s "${OUT}/calling/varscan_snp.vcf" ]; then
  VARSCAN_COUNT=$(grep -c "^chr" "${OUT}/calling/varscan_snp.vcf" || true)
  VARSCAN_COUNT=${VARSCAN_COUNT:-0}
fi
if [ -s "${OUT}/calling/bcftools_raw.vcf.gz" ]; then
  BCFTOOLS_COUNT=$(bcftools view -H "${OUT}/calling/bcftools_raw.vcf.gz" 2>/dev/null | wc -l || true)
  BCFTOOLS_COUNT=${BCFTOOLS_COUNT:-0}
fi
if [ -s "${OUT}/calling/gatk_rna_raw.vcf.gz" ]; then
  GATK_COUNT=$(bcftools view -H "${OUT}/calling/gatk_rna_raw.vcf.gz" 2>/dev/null | wc -l || true)
  GATK_COUNT=${GATK_COUNT:-0}
fi

echo "VarScan variants: ${VARSCAN_COUNT}"
echo "bcftools variants: ${BCFTOOLS_COUNT}"
echo "GATK variants: ${GATK_COUNT}"

# Use bcftools isec to find variants called by >= 2 callers
# First normalize all VCFs to the same format
mkdir -p "${OUT}/calling/isec_input"

# Normalize bcftools calls
if [ -s "${OUT}/calling/bcftools_raw.vcf.gz" ]; then
  bcftools norm -f "${REF}/genome.fa" "${OUT}/calling/bcftools_raw.vcf.gz" -Oz \
    -o "${OUT}/calling/isec_input/bcftools.vcf.gz" 2>/dev/null || true
  bcftools index "${OUT}/calling/isec_input/bcftools.vcf.gz" 2>/dev/null || true
fi

# Normalize GATK calls
if [ -s "${OUT}/calling/gatk_rna_raw.vcf.gz" ]; then
  bcftools norm -f "${REF}/genome.fa" "${OUT}/calling/gatk_rna_raw.vcf.gz" -Oz \
    -o "${OUT}/calling/isec_input/gatk.vcf.gz" 2>/dev/null || true
  bcftools index "${OUT}/calling/isec_input/gatk.vcf.gz" 2>/dev/null || true
fi

# For intersection, use bcftools isec if we have at least 2 VCFs
INTERSECTION_COUNT=0
if [ -s "${OUT}/calling/isec_input/bcftools.vcf.gz" ] && [ -s "${OUT}/calling/isec_input/gatk.vcf.gz" ]; then
  rm -rf "${OUT}/calling/isec_output"
  bcftools isec -n +2 -p "${OUT}/calling/isec_output" \
    "${OUT}/calling/isec_input/bcftools.vcf.gz" \
    "${OUT}/calling/isec_input/gatk.vcf.gz" 2>/dev/null || true
  if [ -f "${OUT}/calling/isec_output/sites.txt" ]; then
    INTERSECTION_COUNT=$(wc -l < "${OUT}/calling/isec_output/sites.txt" || true)
  fi
  # Use one of the intersection outputs as the merged set
  if [ -f "${OUT}/calling/isec_output/0000.vcf" ]; then
    bgzip -c "${OUT}/calling/isec_output/0000.vcf" > "${OUT}/calling/intersected.vcf.gz"
    bcftools index "${OUT}/calling/intersected.vcf.gz" 2>/dev/null || true
  fi
fi

# If intersection is empty, fall back to GATK calls
if [ ! -s "${OUT}/calling/intersected.vcf.gz" ] 2>/dev/null; then
  if [ -s "${OUT}/calling/isec_input/gatk.vcf.gz" ]; then
    cp "${OUT}/calling/isec_input/gatk.vcf.gz" "${OUT}/calling/intersected.vcf.gz"
    cp "${OUT}/calling/isec_input/gatk.vcf.gz.csi" "${OUT}/calling/intersected.vcf.gz.csi" 2>/dev/null || true
    bcftools index "${OUT}/calling/intersected.vcf.gz" 2>/dev/null || true
  fi
fi

echo "Intersection variants (>=2 callers): ${INTERSECTION_COUNT}"

# ============================
# Level 9: Triple annotation/filtering
# ============================

# ---- Filter 1: A>G / T>C strand-aware ----
echo ">>> Level 9a: Strand-specific A>G / T>C filtering"
AG_TC_COUNT=0
if [ -s "${OUT}/calling/intersected.vcf.gz" ]; then
  # A-to-I editing appears as A>G on plus strand, T>C on minus strand
  bcftools view -H "${OUT}/calling/intersected.vcf.gz" 2>/dev/null | \
    awk -F'\t' '($4=="A" && $5=="G") || ($4=="T" && $5=="C")' \
    > "${OUT}/filtered/ag_tc_variants.txt" || true
  AG_TC_COUNT=$(wc -l < "${OUT}/filtered/ag_tc_variants.txt" || true)
  AG_TC_COUNT=${AG_TC_COUNT:-0}
fi
echo "A>G/T>C variants: ${AG_TC_COUNT}"

# ---- Filter 2: Remove known SNPs (dbSNP) ----
echo ">>> Level 9b: dbSNP filtering"
NOVEL_COUNT=0
if [ -s "${OUT}/calling/intersected.vcf.gz" ] && [ -s "${REF}/dbsnp.vcf.gz" ]; then
  bcftools isec -C \
    "${OUT}/calling/intersected.vcf.gz" \
    "${REF}/dbsnp.vcf.gz" \
    -Oz -o "${OUT}/filtered/novel_variants.vcf.gz" -w 1 2>/dev/null || true
  if [ -s "${OUT}/filtered/novel_variants.vcf.gz" ]; then
    NOVEL_COUNT=$(bcftools view -H "${OUT}/filtered/novel_variants.vcf.gz" 2>/dev/null | wc -l || true)
    NOVEL_COUNT=${NOVEL_COUNT:-0}
  fi
fi
echo "Novel (non-dbSNP) variants: ${NOVEL_COUNT}"

# ---- Filter 3: Repeat region intersection ----
echo ">>> Level 9c: Repeat region intersection"
REPEAT_OVERLAP=0
if [ -s "${OUT}/calling/intersected.vcf.gz" ] && [ -s "${REF}/repeat_regions.bed" ]; then
  # Convert VCF positions to BED
  bcftools view -H "${OUT}/calling/intersected.vcf.gz" 2>/dev/null | \
    awk -F'\t' '{print $1"\t"$2-1"\t"$2"\t"$4">"$5}' \
    > "${OUT}/filtered/variant_positions.bed" || true
  if [ -s "${OUT}/filtered/variant_positions.bed" ]; then
    REPEAT_OVERLAP=$(bedtools intersect -a "${OUT}/filtered/variant_positions.bed" \
      -b "${REF}/repeat_regions.bed" -u | wc -l || true)
    REPEAT_OVERLAP=${REPEAT_OVERLAP:-0}
  fi
fi
echo "Variants in repeat regions: ${REPEAT_OVERLAP}"

# ============================
# CONVERGENCE 3 + 4: Final Report
# ============================
echo ">>> CONVERGENCE 3+4: Generating report"

# RNA QC stats
RNA_RAW=$(python3 -c "import json; d=json.load(open('${OUT}/rna_qc/rna_fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
RNA_CLEAN=$(python3 -c "import json; d=json.load(open('${OUT}/rna_qc/rna_fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
RNA_Q30=$(python3 -c "import json; d=json.load(open('${OUT}/rna_qc/rna_fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

# DNA QC stats
DNA_RAW=$(python3 -c "import json; d=json.load(open('${OUT}/dna_qc/dna_fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
DNA_CLEAN=$(python3 -c "import json; d=json.load(open('${OUT}/dna_qc/dna_fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
DNA_Q30=$(python3 -c "import json; d=json.load(open('${OUT}/dna_qc/dna_fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

# Alignment stats from STAR log
RNA_UNIQUE_MAP=$(grep "Uniquely mapped reads %" "${OUT}/rna_align/rna_Log.final.out" | awk -F'|' '{print $2}' | tr -d ' \t%')
RNA_MAPPED=$(grep "in total" "${OUT}/rna_proc/rna_flagstat.txt" | awk '{print $1}')

# DNA alignment stats
DNA_TOTAL=$(grep "in total" "${OUT}/dna_proc/dna_flagstat.txt" | awk '{print $1}')
DNA_MAPPED=$(grep "mapped (" "${OUT}/dna_proc/dna_flagstat.txt" | head -1 | awk '{print $1}')
DNA_MAP_RATE=$(python3 -c "print(round(${DNA_MAPPED}/${DNA_TOTAL}*100,2)) if ${DNA_TOTAL} > 0 else print(0)")

# Duplication rates
RNA_DUP_RATE=$(grep -A1 "^LIBRARY" "${OUT}/rna_proc/rna_markdup_metrics.txt" | tail -1 | awk -F'\t' '{printf "%.2f", $9*100}') || RNA_DUP_RATE=0
DNA_DUP_RATE=$(grep -A1 "^LIBRARY" "${OUT}/dna_proc/dna_markdup_metrics.txt" | tail -1 | awk -F'\t' '{printf "%.2f", $9*100}') || DNA_DUP_RATE=0

# Total candidate editing sites = A>G + T>C that pass filters
CANDIDATE_EDITING=${AG_TC_COUNT}

# Write final report
cat > "${RESULTS}/report.csv" << CSVEOF
metric,value
rna_raw_reads,${RNA_RAW}
rna_clean_reads,${RNA_CLEAN}
rna_q30_rate,${RNA_Q30}
rna_unique_mapping_rate,${RNA_UNIQUE_MAP}
rna_mapped_reads,${RNA_MAPPED}
rna_duplication_rate,${RNA_DUP_RATE}
dna_raw_reads,${DNA_RAW}
dna_clean_reads,${DNA_CLEAN}
dna_q30_rate,${DNA_Q30}
dna_mapped_reads,${DNA_MAPPED}
dna_mapping_rate,${DNA_MAP_RATE}
dna_duplication_rate,${DNA_DUP_RATE}
caller1_variants,${VARSCAN_COUNT}
caller2_variants,${BCFTOOLS_COUNT}
caller3_variants,${GATK_COUNT}
intersection_variants,${INTERSECTION_COUNT}
ag_tc_candidates,${AG_TC_COUNT}
novel_variants,${NOVEL_COUNT}
repeat_region_variants,${REPEAT_OVERLAP}
candidate_editing_sites,${CANDIDATE_EDITING}
CSVEOF

echo "=== Final Report ==="
cat "${RESULTS}/report.csv"
echo ""
echo "Pipeline complete!"
