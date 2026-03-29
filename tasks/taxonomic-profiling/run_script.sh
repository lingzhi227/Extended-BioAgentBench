#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Multi-classifier Taxonomic Profiling Pipeline
# =============================================================================
# DAG (depth=8, convergence=3):
#   fastp(QC) → [classifier1→abundance ‖ classifier2 ‖ classifier3] CONVERGE₁(merge)
#     → [standardize(taxa) ‖ diversity(stats)] CONVERGE₂
#     → consensus_profile CONVERGE₃(multi-tool) → report
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data/fastq"
REFDIR="reference"
OUTDIR="outputs"
RESDIR="results"

mkdir -p "$OUTDIR"/{fastp,kraken2,bracken,centrifuge,kaiju,merged,diversity} "$RESDIR"

SAMPLES=("sample1" "sample2")

# ─── Step 1: Quality control ─────────────────────────────────────────────────
echo "[Step 1] Quality control..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/fastp/${SID}_R1.fastq.gz" ]; then
    fastp \
      -i "$DATADIR/${SID}_R1.fastq.gz" \
      -I "$DATADIR/${SID}_R2.fastq.gz" \
      -o "$OUTDIR/fastp/${SID}_R1.fastq.gz" \
      -O "$OUTDIR/fastp/${SID}_R2.fastq.gz" \
      -j "$OUTDIR/fastp/${SID}_fastp.json" \
      -h "$OUTDIR/fastp/${SID}_fastp.html" \
      -w "$THREADS" --detect_adapter_for_pe \
      --length_required 50 --qualified_quality_phred 20 2>/dev/null
  fi
done

# ─── Step 2a: Classifier 1 — k-mer based (branch 1) ─────────────────────────
echo "[Step 2a] Running classifier 1..."
K2DB="$REFDIR/kraken2/testdb-kraken2"
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/kraken2/${SID}.kreport" ]; then
    kraken2 --db "$K2DB" \
      --paired "$OUTDIR/fastp/${SID}_R1.fastq.gz" "$OUTDIR/fastp/${SID}_R2.fastq.gz" \
      --threads "$THREADS" \
      --output "$OUTDIR/kraken2/${SID}.out" \
      --report "$OUTDIR/kraken2/${SID}.kreport" \
      --quick 2>&1 | tail -3
  fi
done

# Abundance re-estimation
echo "[Step 2a+] Abundance re-estimation..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/bracken/${SID}_S.bracken" ]; then
    bracken -d "$REFDIR/bracken" \
      -i "$OUTDIR/kraken2/${SID}.kreport" \
      -o "$OUTDIR/bracken/${SID}_S.bracken" \
      -w "$OUTDIR/bracken/${SID}_S.breport" \
      -r 150 -l S 2>&1 | tail -2 || true
  fi
done

# ─── Step 2b: Classifier 2 — BWT-based (branch 2) ───────────────────────────
echo "[Step 2b] Running classifier 2..."
CFDB="$REFDIR/centrifuge/test-db-centrifuge/taxprofiler_cf"
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/centrifuge/${SID}.kreport" ]; then
    centrifuge -x "$CFDB" \
      -1 "$OUTDIR/fastp/${SID}_R1.fastq.gz" \
      -2 "$OUTDIR/fastp/${SID}_R2.fastq.gz" \
      -p "$THREADS" \
      -S "$OUTDIR/centrifuge/${SID}.out" \
      --report-file "$OUTDIR/centrifuge/${SID}.report" 2>&1 | tail -3
    centrifuge-kreport -x "$CFDB" "$OUTDIR/centrifuge/${SID}.out" \
      > "$OUTDIR/centrifuge/${SID}.kreport" 2>/dev/null || true
  fi
done

# ─── Step 2c: Classifier 3 — protein-level (branch 3) ───────────────────────
echo "[Step 2c] Running classifier 3..."
KJDB="$REFDIR/kaiju/kaiju"
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/kaiju/${SID}.out" ]; then
    kaiju -z "$THREADS" \
      -t "$KJDB/nodes.dmp" \
      -f "$KJDB/proteins.fmi" \
      -i "$OUTDIR/fastp/${SID}_R1.fastq.gz" \
      -j "$OUTDIR/fastp/${SID}_R2.fastq.gz" \
      -o "$OUTDIR/kaiju/${SID}.out" 2>&1 | tail -3
    kaiju2table -t "$KJDB/nodes.dmp" \
      -n "$KJDB/names.dmp" \
      -r species \
      -o "$OUTDIR/kaiju/${SID}_species.tsv" \
      "$OUTDIR/kaiju/${SID}.out" 2>&1 | tail -2
  fi
done

# ── CONVERGENCE POINTS 1-3: merge, standardize, diversity, report ────────────
echo "[Steps 3-5] Merging, diversity, and report..."
python3 scripts/merge_and_report.py

echo "Pipeline complete. Results in results/report.csv"
