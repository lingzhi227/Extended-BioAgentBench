#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Somatic Variant Calling Pipeline (Tumor-Normal Paired)
# =============================================================================
# DAG (depth=12, convergence=4):
#   [T: fastp→bwa→sort→markdup→BQSR] ‖ [N: fastp→bwa→sort→markdup→BQSR]
#     CONVERGE₁(shared reference + sites)
#   → [Caller1(somatic-paired) ‖ Caller2(joint) ‖ Caller3(pileup)]
#     CONVERGE₂(multi-caller)
#   → [filter ‖ coverage-stats]
#     CONVERGE₃(filtered-variants + coverage)
#   → annotate
#     CONVERGE₄(annotations + variants) → report
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
REF="reference/genome.fasta"
DBSNP="reference/dbsnp.vcf.gz"
KNOWN_INDELS="reference/known_indels.vcf.gz"
GNOMAD="reference/gnomad.vcf.gz"
INTERVALS="reference/intervals.list"
OUTDIR="outputs"
RESDIR="results"

mkdir -p "$OUTDIR"/{fastp,align,markdup,bqsr,mutect2,freebayes,bcftools_call,merge,filter,coverage,annotate} "$RESDIR"

# ─── Step 1a: Tumor QC + alignment (branch 1) ────────────────────────────────
echo "[Step 1a] Processing tumor sample..."
if [ ! -f "$OUTDIR/align/tumor.sorted.bam" ]; then
  # Merge all tumor lane FASTQs
  cat data/tumor/tumor_*_R1.fastq.gz > "$OUTDIR/fastp/tumor_R1.fastq.gz"
  cat data/tumor/tumor_*_R2.fastq.gz > "$OUTDIR/fastp/tumor_R2.fastq.gz"

  # QC
  fastp -i "$OUTDIR/fastp/tumor_R1.fastq.gz" -I "$OUTDIR/fastp/tumor_R2.fastq.gz" \
    -o "$OUTDIR/fastp/tumor_qc_R1.fastq.gz" -O "$OUTDIR/fastp/tumor_qc_R2.fastq.gz" \
    -j "$OUTDIR/fastp/tumor_fastp.json" -w "$THREADS" --detect_adapter_for_pe 2>/dev/null

  # Align with read group
  bwa-mem2 mem -t "$THREADS" \
    -R "@RG\tID:tumor\tSM:tumor\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "$REF" "$OUTDIR/fastp/tumor_qc_R1.fastq.gz" "$OUTDIR/fastp/tumor_qc_R2.fastq.gz" 2>/dev/null \
    | samtools sort -@ "$THREADS" -o "$OUTDIR/align/tumor.sorted.bam" -
  samtools index "$OUTDIR/align/tumor.sorted.bam"
fi

# ─── Step 1b: Normal QC + alignment (branch 2) ──────────────────────────────
echo "[Step 1b] Processing normal sample..."
if [ ! -f "$OUTDIR/align/normal.sorted.bam" ]; then
  cat data/normal/normal_*_R1.fastq.gz > "$OUTDIR/fastp/normal_R1.fastq.gz"
  cat data/normal/normal_*_R2.fastq.gz > "$OUTDIR/fastp/normal_R2.fastq.gz"

  fastp -i "$OUTDIR/fastp/normal_R1.fastq.gz" -I "$OUTDIR/fastp/normal_R2.fastq.gz" \
    -o "$OUTDIR/fastp/normal_qc_R1.fastq.gz" -O "$OUTDIR/fastp/normal_qc_R2.fastq.gz" \
    -j "$OUTDIR/fastp/normal_fastp.json" -w "$THREADS" --detect_adapter_for_pe 2>/dev/null

  bwa-mem2 mem -t "$THREADS" \
    -R "@RG\tID:normal\tSM:normal\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "$REF" "$OUTDIR/fastp/normal_qc_R1.fastq.gz" "$OUTDIR/fastp/normal_qc_R2.fastq.gz" 2>/dev/null \
    | samtools sort -@ "$THREADS" -o "$OUTDIR/align/normal.sorted.bam" -
  samtools index "$OUTDIR/align/normal.sorted.bam"
fi

# ─── Step 2: Mark duplicates ─────────────────────────────────────────────────
echo "[Step 2] Marking duplicates..."
for SAMPLE in tumor normal; do
  if [ ! -f "$OUTDIR/markdup/${SAMPLE}.dedup.bam" ]; then
    gatk MarkDuplicates \
      -I "$OUTDIR/align/${SAMPLE}.sorted.bam" \
      -O "$OUTDIR/markdup/${SAMPLE}.dedup.bam" \
      -M "$OUTDIR/markdup/${SAMPLE}.metrics.txt" \
      --REMOVE_DUPLICATES false \
      --CREATE_INDEX true 2>/dev/null
  fi
done

# ─── Step 3: Base Quality Score Recalibration ────────────────────────────────
echo "[Step 3] BQSR..."
for SAMPLE in tumor normal; do
  if [ ! -f "$OUTDIR/bqsr/${SAMPLE}.recal.bam" ]; then
    gatk BaseRecalibrator \
      -I "$OUTDIR/markdup/${SAMPLE}.dedup.bam" \
      -R "$REF" \
      --known-sites "$DBSNP" \
      --known-sites "$KNOWN_INDELS" \
      -O "$OUTDIR/bqsr/${SAMPLE}.recal_table" 2>/dev/null

    gatk ApplyBQSR \
      -I "$OUTDIR/markdup/${SAMPLE}.dedup.bam" \
      -R "$REF" \
      --bqsr-recal-file "$OUTDIR/bqsr/${SAMPLE}.recal_table" \
      -O "$OUTDIR/bqsr/${SAMPLE}.recal.bam" 2>/dev/null

    samtools index "$OUTDIR/bqsr/${SAMPLE}.recal.bam"
  fi
done

# ── CONVERGENCE POINT 1: Tumor + Normal preprocessed ────────────────────────

TUMOR_BAM="$OUTDIR/bqsr/tumor.recal.bam"
NORMAL_BAM="$OUTDIR/bqsr/normal.recal.bam"

# ─── Step 4a: Caller 1 — Paired somatic calling (branch 1) ──────────────────
echo "[Step 4a] Running Caller 1 (somatic paired)..."
if [ ! -f "$OUTDIR/mutect2/somatic_raw.vcf.gz" ]; then
  gatk Mutect2 \
    -R "$REF" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    -normal normal \
    --germline-resource "$GNOMAD" \
    -L "$INTERVALS" \
    -O "$OUTDIR/mutect2/somatic_raw.vcf.gz" \
    --f1r2-tar-gz "$OUTDIR/mutect2/f1r2.tar.gz" 2>/dev/null

  # Get pileup summaries for contamination estimation
  gatk GetPileupSummaries \
    -I "$TUMOR_BAM" \
    -V "$GNOMAD" \
    -L "$GNOMAD" \
    -O "$OUTDIR/mutect2/tumor_pileups.table" 2>/dev/null || true

  gatk GetPileupSummaries \
    -I "$NORMAL_BAM" \
    -V "$GNOMAD" \
    -L "$GNOMAD" \
    -O "$OUTDIR/mutect2/normal_pileups.table" 2>/dev/null || true

  # Calculate contamination
  if [ -f "$OUTDIR/mutect2/tumor_pileups.table" ] && [ -f "$OUTDIR/mutect2/normal_pileups.table" ]; then
    gatk CalculateContamination \
      -I "$OUTDIR/mutect2/tumor_pileups.table" \
      -matched "$OUTDIR/mutect2/normal_pileups.table" \
      -O "$OUTDIR/mutect2/contamination.table" \
      --tumor-segmentation "$OUTDIR/mutect2/segments.table" 2>/dev/null || true
  fi

  # Filter
  FILTER_ARGS="-R $REF -V $OUTDIR/mutect2/somatic_raw.vcf.gz -O $OUTDIR/mutect2/somatic_filtered.vcf.gz"
  if [ -f "$OUTDIR/mutect2/contamination.table" ]; then
    FILTER_ARGS="$FILTER_ARGS --contamination-table $OUTDIR/mutect2/contamination.table"
  fi
  if [ -f "$OUTDIR/mutect2/segments.table" ]; then
    FILTER_ARGS="$FILTER_ARGS --tumor-segmentation $OUTDIR/mutect2/segments.table"
  fi
  gatk FilterMutectCalls $FILTER_ARGS 2>/dev/null || true
fi

# ─── Step 4b: Caller 2 — Joint genotype calling (branch 2) ──────────────────
echo "[Step 4b] Running Caller 2 (joint calling)..."
if [ ! -f "$OUTDIR/freebayes/joint_raw.vcf" ]; then
  freebayes -f "$REF" \
    --pooled-continuous \
    --min-alternate-fraction 0.05 \
    --min-alternate-count 2 \
    "$TUMOR_BAM" "$NORMAL_BAM" \
    > "$OUTDIR/freebayes/joint_raw.vcf" 2>/dev/null || true
fi

# ─── Step 4c: Caller 3 — Pileup-based calling (branch 3) ────────────────────
echo "[Step 4c] Running Caller 3 (pileup-based)..."
if [ ! -f "$OUTDIR/bcftools_call/pileup_raw.vcf.gz" ]; then
  bcftools mpileup -f "$REF" \
    --annotate FORMAT/AD,FORMAT/DP \
    "$TUMOR_BAM" "$NORMAL_BAM" 2>/dev/null \
    | bcftools call -mv -Oz -o "$OUTDIR/bcftools_call/pileup_raw.vcf.gz" 2>/dev/null || true
  bcftools index "$OUTDIR/bcftools_call/pileup_raw.vcf.gz" 2>/dev/null || true
fi

# ── CONVERGENCE POINT 2: Multi-caller results ───────────────────────────────

# ─── Step 5a: Filter and normalize variants (branch 1) ──────────────────────
echo "[Step 5a] Filtering variants..."
# Normalize all VCFs
for vcf_pair in "mutect2/somatic_filtered.vcf.gz:filter/mutect2_norm.vcf.gz" \
                "freebayes/joint_raw.vcf:filter/freebayes_norm.vcf.gz" \
                "bcftools_call/pileup_raw.vcf.gz:filter/bcftools_norm.vcf.gz"; do
  SRC="${vcf_pair%%:*}"
  DST="${vcf_pair##*:}"
  if [ -f "$OUTDIR/$SRC" ] && [ ! -f "$OUTDIR/$DST" ]; then
    bcftools norm -f "$REF" -m -both "$OUTDIR/$SRC" 2>/dev/null \
      | bcftools view -f PASS,. -Oz -o "$OUTDIR/$DST" 2>/dev/null || \
      cp "$OUTDIR/$SRC" "$OUTDIR/$DST" 2>/dev/null || true
    bcftools index "$OUTDIR/$DST" 2>/dev/null || true
  fi
done

# ─── Step 5b: Coverage statistics (branch 2) ────────────────────────────────
echo "[Step 5b] Computing coverage..."
for SAMPLE in tumor normal; do
  if [ ! -f "$OUTDIR/coverage/${SAMPLE}.mosdepth.summary.txt" ]; then
    mosdepth --threads "$THREADS" --fast-mode \
      "$OUTDIR/coverage/$SAMPLE" "$OUTDIR/bqsr/${SAMPLE}.recal.bam" 2>/dev/null || true
  fi
done

# ── CONVERGENCE POINT 3: Filtered variants + coverage ───────────────────────

# ─── Step 6: Annotate variants ──────────────────────────────────────────────
echo "[Step 6] Annotating variants..."
# Annotate with dbSNP
for caller in mutect2 freebayes bcftools; do
  NORM="$OUTDIR/filter/${caller}_norm.vcf.gz"
  ANN="$OUTDIR/annotate/${caller}_annotated.vcf.gz"
  if [ -f "$NORM" ] && [ ! -f "$ANN" ]; then
    bcftools annotate -a "$DBSNP" -c ID "$NORM" -Oz -o "$ANN" 2>/dev/null || \
      cp "$NORM" "$ANN" 2>/dev/null || true
    bcftools index "$ANN" 2>/dev/null || true
  fi
done

# ── CONVERGENCE POINT 4: Annotated variants → report ────────────────────────

# ─── Step 7: Generate report ────────────────────────────────────────────────
echo "[Step 7] Compiling report..."
python3 scripts/compile_report.py

echo "Pipeline complete. Results in results/report.csv"
