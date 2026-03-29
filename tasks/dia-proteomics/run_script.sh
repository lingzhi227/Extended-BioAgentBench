#!/usr/bin/env bash
set -euo pipefail
#
# DIA Proteomics Quantification Pipeline
#
# DAG Structure (depth=9, convergence=2):
#
#   mzML_1 → openswath(1a) ─┐
#   mzML_2 → openswath(1b) ─┤
#   mzML_3 → openswath(1c) ─┤
#   mzML_4 → openswath(1d) ─┤
#                            ↓ CONVERGE1 (merge per-run results)
#                     pyprophet_merge(2)
#                       → pyprophet_score(3)
#                         → pyprophet_export(4)
#                           ├→ feature_alignment(5a) ──┐
#                           └→ peptide_summary(5b) ────┤ CONVERGE2
#                                                       ↓
#                                               msstats(6) → report(7)
#
# Inputs:  data/*.mzML, reference/spectral_library.pqp, reference/irts.pqp,
#          reference/swath_windows.txt, reference/sample_sheet.tsv
# Outputs: results/report.csv

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
REF_DIR="${SCRIPT_DIR}/reference"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

mkdir -p "${OUTPUT_DIR}"/{openswath,pyprophet,tric,msstats} "${RESULTS_DIR}"

echo "=== DIA Proteomics Quantification Pipeline ==="
echo "Threads: ${THREADS}"

# ============================================================
# Step 1: OpenSwathWorkflow (per-run DIA analysis)
# ============================================================
MZML_FILES=$(ls "${DATA_DIR}"/DIA_test_*.mzML 2>/dev/null)
N_RUNS=$(echo "$MZML_FILES" | wc -l)
echo "Found ${N_RUNS} DIA runs"

for mzml in ${MZML_FILES}; do
    BASENAME=$(basename "$mzml" .mzML)
    TSV_OUT="${OUTPUT_DIR}/openswath/${BASENAME}.tsv"
    if [ ! -f "${TSV_OUT}" ]; then
        echo "[Step 1] OpenSwath: ${BASENAME}..."
        OpenSwathWorkflow \
            -in "$mzml" \
            -tr "${REF_DIR}/spectral_library.pqp" \
            -tr_irt "${REF_DIR}/irts.pqp" \
            -out_tsv "${TSV_OUT}" \
            -threads ${THREADS} \
            -sort_swath_maps \
            -force \
            -readOptions cacheWorkingInMemory \
            -batchSize 250 \
            -min_upper_edge_dist 1 \
            -Scoring:stop_report_after_feature 5 \
            -rt_extraction_window 600 \
            -mz_extraction_window 30 \
            -mz_extraction_window_unit ppm \
            -RTNormalization:estimateBestPeptides \
            -RTNormalization:outlierMethod none \
            -RTNormalization:NrRTBins 2 \
            -RTNormalization:MinBinsFilled 1 \
            -RTNormalization:MinPeptidesPerBin 1 \
            -RTNormalization:InitialQualityCutoff -2 \
            -RTNormalization:OverallQualityCutoff 0 2>&1 | tail -5
    else
        echo "[Step 1] ${BASENAME} already processed, skipping."
    fi
done

# ============================================================
# Step 2-3: Score each run with pyprophet (CONVERGE1 = all runs scored)
# ============================================================
for tsv in "${OUTPUT_DIR}"/openswath/*.tsv; do
    BASENAME=$(basename "$tsv" .tsv)
    SCORED="${OUTPUT_DIR}/pyprophet/${BASENAME}_scored.tsv"
    if [ ! -f "${SCORED}" ]; then
        echo "[Step 2-3] Scoring: ${BASENAME}..."
        cp "$tsv" "${OUTPUT_DIR}/pyprophet/${BASENAME}.tsv"
        cd "${OUTPUT_DIR}/pyprophet"
        pyprophet --target.overwrite --ignore.invalid_score_columns \
            "${BASENAME}.tsv" 2>&1 | tail -5
        # pyprophet creates _with_dscore_filtered.csv
        FILTERED=$(ls "${BASENAME}"*_with_dscore_filtered.csv 2>/dev/null | head -1)
        if [ -n "$FILTERED" ] && [ -f "$FILTERED" ]; then
            mv "$FILTERED" "${SCORED}"
        else
            DSCORE=$(ls "${BASENAME}"*_with_dscore.csv 2>/dev/null | head -1)
            if [ -n "$DSCORE" ] && [ -f "$DSCORE" ]; then
                mv "$DSCORE" "${SCORED}"
            fi
        fi
        cd "${SCRIPT_DIR}"
    else
        echo "[Step 2-3] ${BASENAME} already scored, skipping."
    fi
done

# ============================================================
# Step 4: Combine scored results for alignment
# ============================================================
EXPORT_TSV="${OUTPUT_DIR}/pyprophet/combined.tsv"
if [ ! -f "${EXPORT_TSV}" ]; then
    echo "[Step 4] Combining scored results..."
    # Concatenate all scored TSVs (keep header from first, skip headers from rest)
    FIRST=true
    for scored in "${OUTPUT_DIR}"/pyprophet/*_scored.tsv; do
        if [ "$FIRST" = true ]; then
            cat "$scored" > "${EXPORT_TSV}"
            FIRST=false
        else
            tail -n +2 "$scored" >> "${EXPORT_TSV}"
        fi
    done
    echo "  Combined: $(wc -l < "${EXPORT_TSV}") lines"
else
    echo "[Step 4] Combining already done, skipping."
fi

# ============================================================
# Step 5a: Feature alignment across runs (TRIC)
# ============================================================
ALIGNED="${OUTPUT_DIR}/tric/aligned.tsv"
if [ ! -f "${ALIGNED}" ]; then
    echo "[Step 5a] Aligning features across runs (TRIC)..."
    feature_alignment.py \
        --in "${EXPORT_TSV}" \
        --out "${ALIGNED}" \
        --method best_overall \
        --max_rt_diff 300 \
        --alignment_score 0.01 \
        --fdr_cutoff 0.05 2>&1 | tail -10
else
    echo "[Step 5a] Alignment already done, skipping."
fi

# ============================================================
# Step 5b + 6: MSstats quantification (CONVERGE2)
# ============================================================
MSSTATS_OUT="${OUTPUT_DIR}/msstats"
if [ ! -f "${MSSTATS_OUT}/msstats_results.csv" ]; then
    echo "[Step 6] Running MSstats..."
    Rscript "${SCRIPT_DIR}/scripts/msstats_analysis.R" \
        "${ALIGNED}" "${REF_DIR}/sample_sheet.tsv" "${MSSTATS_OUT}"
else
    echo "[Step 6] MSstats already done, skipping."
fi

# ============================================================
# Step 7: Generate report
# ============================================================
echo "[Step 7] Generating report..."

# Count peptides and proteins from export
if [ -f "${EXPORT_TSV}" ]; then
    N_PEPTIDES=$(awk -F'\t' 'NR>1 {print $0}' "${EXPORT_TSV}" | cut -f1 | sort -u | wc -l || echo "0")
    N_PROTEINS=$(awk -F'\t' 'NR>1' "${EXPORT_TSV}" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /ProteinName/) {col=i; break}} NR>1 {print $col}' "${EXPORT_TSV}" | sort -u | wc -l || echo "0")
    # Better approach: check header
    HEADER=$(head -1 "${EXPORT_TSV}")
    PROT_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n "ProteinName" | cut -d: -f1 || echo "0")
    PEP_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n "^Sequence$\|FullPeptideName\|PeptideSequence" | head -1 | cut -d: -f1 || echo "0")
    if [ "$PROT_COL" -gt 0 ]; then
        N_PROTEINS=$(awk -F'\t' -v c=$PROT_COL 'NR>1 {print $c}' "${EXPORT_TSV}" | sort -u | wc -l)
    fi
    if [ "$PEP_COL" -gt 0 ]; then
        N_PEPTIDES=$(awk -F'\t' -v c=$PEP_COL 'NR>1 {print $c}' "${EXPORT_TSV}" | sort -u | wc -l)
    fi
    N_TOTAL_FEATURES=$(tail -n +2 "${EXPORT_TSV}" | wc -l)
else
    N_PEPTIDES=0
    N_PROTEINS=0
    N_TOTAL_FEATURES=0
fi

# Count aligned features
N_ALIGNED=0
if [ -f "${ALIGNED}" ]; then
    N_ALIGNED=$(tail -n +2 "${ALIGNED}" | wc -l)
fi

# MSstats results
N_DE_PROTEINS=0
if [ -f "${MSSTATS_OUT}/msstats_results.csv" ]; then
    N_DE_PROTEINS=$(awk -F',' 'NR>1 && $NF < 0.05 {count++} END {print count+0}' "${MSSTATS_OUT}/msstats_results.csv" || echo "0")
fi

cat > "${RESULTS_DIR}/report.csv" << CSVEOF
metric,value
dia_runs_processed,${N_RUNS}
total_features,${N_TOTAL_FEATURES}
peptides_identified,${N_PEPTIDES}
proteins_identified,${N_PROTEINS}
aligned_features,${N_ALIGNED}
differentially_expressed_proteins,${N_DE_PROTEINS}
fdr_threshold,0.01
CSVEOF

echo ""
echo "=== Final Report ==="
cat "${RESULTS_DIR}/report.csv"

N_EMPTY=$(grep -cE ",,|,$|,NA$" "${RESULTS_DIR}/report.csv" || true)
echo ""
echo "Empty/NA values: ${N_EMPTY}"
echo "=== Pipeline complete ==="
