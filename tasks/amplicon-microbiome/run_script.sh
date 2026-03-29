#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 16S Amplicon Microbiome Analysis Pipeline
# =============================================================================
# DAG (depth=9, convergence=3):
#   cutadapt(primers) → filter → learnErrors → denoise → mergePairs → removeChimeras
#     → [assignTaxonomy ‖ MAFFT→FastTree] CONVERGE₁(taxonomy+phylogeny)
#     → [functionalPrediction ‖ diversityMetrics] CONVERGE₂(function+diversity)
#     → differentialAbundance CONVERGE₃(metadata+results) → report
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data"
REFDIR="reference"
OUTDIR="outputs"
RESDIR="results"

mkdir -p "$OUTDIR"/{trimmed,filtered,dada2,taxonomy,phylogeny,picrust2,diversity} "$RESDIR"

# Primer sequences (515F / 806R, V4 region)
FWD_PRIMER="GTGYCAGCMGCCGCGGTAA"
REV_PRIMER="GGACTACNVGGGTWTCTAAT"
REV_PRIMER_RC="ATTAGAWACCCBDGTAGTCC"
FWD_PRIMER_RC="TTACCGCGGCKGCTGRCAC"

SAMPLES=("sample1:1_S103" "sample1a:1a_S103" "sample2:2_S115" "sample2a:2a_S115")

# ─── Step 1: Primer removal with cutadapt ────────────────────────────────────
echo "[Step 1] Removing primers..."
for entry in "${SAMPLES[@]}"; do
  SID="${entry%%:*}"
  PREFIX="${entry##*:}"
  R1="$DATADIR/${PREFIX}_L001_R1_001.fastq.gz"
  R2="$DATADIR/${PREFIX}_L001_R2_001.fastq.gz"
  O1="$OUTDIR/trimmed/${SID}_R1.fastq.gz"
  O2="$OUTDIR/trimmed/${SID}_R2.fastq.gz"
  if [ ! -f "$O1" ]; then
    cutadapt \
      -g "$FWD_PRIMER" -G "$REV_PRIMER" \
      -a "$REV_PRIMER_RC" -A "$FWD_PRIMER_RC" \
      --discard-untrimmed --minimum-length 50 \
      -j "$THREADS" \
      -o "$O1" -p "$O2" "$R1" "$R2" \
      > "$OUTDIR/trimmed/${SID}_cutadapt.log" 2>&1
  fi
done

# ─── Steps 2-6: DADA2 denoising pipeline ─────────────────────────────────────
echo "[Steps 2-6] Running denoising pipeline..."
if [ ! -f "$OUTDIR/dada2/asv_table.tsv" ]; then
  Rscript --vanilla scripts/01_dada2.R
fi

# ─── Step 7a: Taxonomy classification (branch 1) ─────────────────────────────
echo "[Step 7a] Classifying taxonomy..."
if [ ! -f "$OUTDIR/taxonomy/taxonomy.tsv" ]; then
  Rscript --vanilla scripts/02_taxonomy.R
fi

# ─── Step 7b: Phylogenetic tree (branch 2) ───────────────────────────────────
echo "[Step 7b] Building phylogenetic tree..."
if [ ! -f "$OUTDIR/phylogeny/tree.nwk" ]; then
  mafft --auto --thread "$THREADS" \
    "$OUTDIR/dada2/asv_seqs.fasta" > "$OUTDIR/phylogeny/aligned.fasta" 2>/dev/null
  FastTree -gtr -nt -fastest \
    "$OUTDIR/phylogeny/aligned.fasta" > "$OUTDIR/phylogeny/tree.nwk" 2>/dev/null
fi

# ── CONVERGENCE POINT 1: taxonomy + phylogeny ────────────────────────────────

# ─── Step 8a: Functional prediction (branch 1) ───────────────────────────────
echo "[Step 8a] Predicting functional potential..."
if [ ! -d "$OUTDIR/picrust2/out/pathways_out" ]; then
  # Prepare frequency table
  python3 << 'PYEOF'
import csv
with open("outputs/dada2/asv_table.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    samples, data = [], {}
    for row in reader:
        sid = row["sample_id"]
        samples.append(sid)
        for k, v in row.items():
            if k != "sample_id":
                data.setdefault(k, {})[sid] = int(v)
asvs = sorted(data.keys(), key=lambda x: int(x.replace("ASV","")))
with open("outputs/picrust2/freq_table.tsv", "w") as f:
    f.write("#OTU ID\t" + "\t".join(samples) + "\n")
    for asv in asvs:
        vals = [str(data[asv].get(s, 0)) for s in samples]
        f.write(asv + "\t" + "\t".join(vals) + "\n")
PYEOF

  rm -rf "$OUTDIR/picrust2/out"
  picrust2_pipeline.py \
    -s "$OUTDIR/dada2/asv_seqs.fasta" \
    -i "$OUTDIR/picrust2/freq_table.tsv" \
    -o "$OUTDIR/picrust2/out" \
    -p "$THREADS" --verbose 2>&1 | tail -5 || true
fi

# ─── Step 8b: Diversity metrics (branch 2) ───────────────────────────────────
echo "[Step 8b] Calculating diversity metrics..."
if [ ! -f "$OUTDIR/diversity/diversity_metrics.tsv" ]; then
  Rscript --vanilla scripts/03_diversity.R
fi

# ── CONVERGENCE POINT 2: functional prediction + diversity ───────────────────

# ─── Step 9: Differential abundance testing ──────────────────────────────────
echo "[Step 9] Testing differential abundance..."
if [ ! -f "$OUTDIR/diversity/da_results.tsv" ]; then
  Rscript --vanilla scripts/04_da.R
fi

# ── CONVERGENCE POINT 3: DA + metadata + all results ────────────────────────

# ─── Step 10: Compile final report ───────────────────────────────────────────
echo "[Step 10] Compiling report..."
Rscript --vanilla scripts/05_report.R

echo "Pipeline complete. Results in results/report.csv"
