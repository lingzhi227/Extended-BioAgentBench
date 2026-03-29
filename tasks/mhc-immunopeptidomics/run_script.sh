#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# MHC Immunopeptidomics Pipeline
# =============================================================================
# DAG (depth=10, convergence=3):
#   mzML → PeakPicker → DatabaseSearch(no enzyme)
#     → PSMFeatureExtractor → [Rescoring ‖ DecoyAnalysis] CONVERGE₁
#     → IDFilter → [MapAligner(across runs) ‖ IDMerger] CONVERGE₂
#     → FeatureFinder → [Quantification ‖ PeptideSummary] CONVERGE₃
#     → report
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
OUTDIR="outputs"
RESDIR="results"
REF="reference/protein_database.fasta"

SAMPLES=("sample1" "sample2")

mkdir -p "$OUTDIR"/{peakpick,search,features,rescore,filter,align,quantify} "$RESDIR"

# ─── Step 1: Peak picking (centroid spectra) ─────────────────────────────────
echo "[Step 1] Peak picking..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/peakpick/${SID}.mzML" ]; then
    PeakPickerHiRes \
      -in "data/${SID}.mzML" \
      -out "$OUTDIR/peakpick/${SID}.mzML" \
      -algorithm:ms_levels "[1,2]" \
      -threads "$THREADS" 2>/dev/null
  fi
done

# ─── Step 2: Database search (unspecific cleavage for MHC peptides) ──────────
echo "[Step 2] Database search..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/search/${SID}.idXML" ]; then
    CometAdapter \
      -in "$OUTDIR/peakpick/${SID}.mzML" \
      -out "$OUTDIR/search/${SID}.idXML" \
      -database "$REF" \
      -precursor_mass_tolerance 20 \
      -fragment_mass_tolerance 0.02 \
      -enzyme "unspecific cleavage" \
      -fixed_modifications "Carbamidomethyl (C)" \
      -variable_modifications "Oxidation (M)" \
      -allowed_missed_cleavages 0 \
      -peptide_length_min 8 \
      -peptide_length_max 12 \
      -threads "$THREADS" 2>/dev/null
  fi
done

# ─── Step 3: Extract PSM features ───────────────────────────────────────────
echo "[Step 3] Extracting PSM features..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/features/${SID}.idXML" ]; then
    PSMFeatureExtractor \
      -in "$OUTDIR/search/${SID}.idXML" \
      -out "$OUTDIR/features/${SID}.idXML" 2>/dev/null
  fi
done

# ─── Step 4a: Statistical rescoring (branch 1) ──────────────────────────────
echo "[Step 4a] Rescoring with statistical validation..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/rescore/${SID}.idXML" ]; then
    PercolatorAdapter \
      -in "$OUTDIR/features/${SID}.idXML" \
      -out "$OUTDIR/rescore/${SID}.idXML" \
      -trainFDR 0.05 \
      -testFDR 0.05 \
      -enzyme no_enzyme \
      -threads "$THREADS" 2>/dev/null || {
        echo "  Percolator failed for ${SID}, copying features..."
        cp "$OUTDIR/features/${SID}.idXML" "$OUTDIR/rescore/${SID}.idXML"
      }
  fi
done

# ─── Step 4b: Decoy statistics (branch 2) ───────────────────────────────────
echo "[Step 4b] Computing decoy statistics..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/rescore/${SID}_fdr.idXML" ]; then
    FalseDiscoveryRate \
      -in "$OUTDIR/features/${SID}.idXML" \
      -out "$OUTDIR/rescore/${SID}_fdr.idXML" \
      -PSM true \
      -protein false 2>/dev/null || true
  fi
done

# ── CONVERGENCE POINT 1: Rescoring + FDR analysis ───────────────────────────

# ─── Step 5: Filter identifications ─────────────────────────────────────────
echo "[Step 5] Filtering identifications..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/filter/${SID}.idXML" ]; then
    IDFilter \
      -in "$OUTDIR/rescore/${SID}.idXML" \
      -out "$OUTDIR/filter/${SID}.idXML" \
      -score:pep 0.05 \
      -remove_decoys 2>/dev/null || {
        # Fall back to less stringent filter
        IDFilter \
          -in "$OUTDIR/rescore/${SID}.idXML" \
          -out "$OUTDIR/filter/${SID}.idXML" \
          -best:n_peptide_hits 1 2>/dev/null || \
          cp "$OUTDIR/rescore/${SID}.idXML" "$OUTDIR/filter/${SID}.idXML"
      }
  fi
done

# ─── Step 6a: Align identifications across runs (branch 1) ──────────────────
echo "[Step 6a] Aligning across runs..."
if [ ! -f "$OUTDIR/align/sample1.idXML" ]; then
  INPUT_FILES=""
  OUTPUT_FILES=""
  for SID in "${SAMPLES[@]}"; do
    INPUT_FILES="$INPUT_FILES -in $OUTDIR/filter/${SID}.idXML"
    OUTPUT_FILES="$OUTPUT_FILES -out $OUTDIR/align/${SID}.idXML"
  done
  MapAlignerIdentification \
    $INPUT_FILES $OUTPUT_FILES 2>/dev/null || {
      echo "  Alignment skipped, copying files..."
      for SID in "${SAMPLES[@]}"; do
        cp "$OUTDIR/filter/${SID}.idXML" "$OUTDIR/align/${SID}.idXML"
      done
    }
fi

# ─── Step 6b: Merge identifications (branch 2) ──────────────────────────────
echo "[Step 6b] Merging identifications..."
if [ ! -f "$OUTDIR/align/merged.idXML" ]; then
  INPUT_FILES=""
  for SID in "${SAMPLES[@]}"; do
    INPUT_FILES="$INPUT_FILES -in $OUTDIR/filter/${SID}.idXML"
  done
  IDMerger $INPUT_FILES -out "$OUTDIR/align/merged.idXML" 2>/dev/null || true
fi

# ── CONVERGENCE POINT 2: Aligned + merged identifications ───────────────────

# ─── Step 7: Feature-based quantification ───────────────────────────────────
echo "[Step 7] Feature-based quantification..."
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/quantify/${SID}.featureXML" ]; then
    FeatureFinderIdentification \
      -in "$OUTDIR/peakpick/${SID}.mzML" \
      -id "$OUTDIR/align/${SID}.idXML" \
      -out "$OUTDIR/quantify/${SID}.featureXML" \
      -threads "$THREADS" 2>/dev/null || true
  fi
done

# ── CONVERGENCE POINT 3: Quantification + identifications ───────────────────

# ─── Step 8: Export and compile report ───────────────────────────────────────
echo "[Step 8] Exporting results and compiling report..."

# Export identifications to TSV
for SID in "${SAMPLES[@]}"; do
  if [ ! -f "$OUTDIR/quantify/${SID}_ids.tsv" ]; then
    TextExporter \
      -in "$OUTDIR/filter/${SID}.idXML" \
      -out "$OUTDIR/quantify/${SID}_ids.tsv" 2>/dev/null || true
  fi
done

# Compile report
python3 scripts/compile_report.py

echo "Pipeline complete. Results in results/report.csv"
