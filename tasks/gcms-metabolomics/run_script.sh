#!/usr/bin/env bash
set -euo pipefail
#
# GC-MS Metabolomics Profiling Pipeline
#
# DAG Structure (depth=9, convergence=2):
#
#   read_mzml(1) → peak_detect(2) → group_peaks(3) → align_rt(4) → regroup(5) → fill_peaks(6)
#     ├─→ annotate_groups(7a) ──────────────────┐
#     └─→ compute_retention_index(7b) ──────────┤ CONVERGE1
#                                                ↓
#     ├─→ spectral_library_match(8a) ───────────┐
#     └─→ statistical_summary(8b) ──────────────┤ CONVERGE2
#                                                ↓
#                                          report(9)
#
# Inputs:  data/*.mzML, data/spectral_library.msp, data/reference_alkanes.csv, data/sampleMetadata.tsv
# Outputs: results/report.csv

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

mkdir -p "${OUTPUT_DIR}"/{xcms,camera,ri,matching} "${RESULTS_DIR}"

echo "=== GC-MS Metabolomics Pipeline ==="
echo "Threads: ${THREADS}"
echo "Data: ${DATA_DIR}"

# ============================================================
# Steps 1-6: XCMS processing (read → peaks → group → align → regroup → fill)
# ============================================================
if [ ! -f "${OUTPUT_DIR}/xcms/xdata.rds" ]; then
    echo "[Step 1-6] Running XCMS processing..."
    Rscript "${SCRIPT_DIR}/scripts/xcms_process.R" \
        "${DATA_DIR}" "${OUTPUT_DIR}/xcms" "${THREADS}"
else
    echo "[Step 1-6] XCMS processing already done, skipping."
fi

# Verify XCMS output
if [ ! -f "${OUTPUT_DIR}/xcms/peak_matrix.csv" ] || [ ! -f "${OUTPUT_DIR}/xcms/feature_defs.csv" ]; then
    echo "ERROR: XCMS processing failed - missing output files"
    exit 1
fi

# ============================================================
# Step 7a: CAMERA annotation (isotopes + adducts + pseudospectra)
# ============================================================
if [ ! -f "${OUTPUT_DIR}/camera/camera_summary.csv" ]; then
    echo "[Step 7a] Running spectral group annotation..."
    Rscript "${SCRIPT_DIR}/scripts/camera_annotate.R" \
        "${OUTPUT_DIR}/xcms" "${OUTPUT_DIR}/camera"
else
    echo "[Step 7a] Annotation already done, skipping."
fi

# ============================================================
# Step 7b: Retention index calculation (parallel with 7a)
# ============================================================
if [ ! -f "${OUTPUT_DIR}/ri/ri_summary.csv" ]; then
    echo "[Step 7b] Computing retention indices..."
    python3 "${SCRIPT_DIR}/scripts/ri_assign.py" \
        "${OUTPUT_DIR}/xcms/feature_defs.csv" \
        "${DATA_DIR}/reference_alkanes.csv" \
        "${OUTPUT_DIR}/ri"
else
    echo "[Step 7b] RI calculation already done, skipping."
fi

# --- CONVERGE1: both annotation and RI complete ---

# ============================================================
# Step 8a: Spectral library matching
# ============================================================
if [ ! -f "${OUTPUT_DIR}/matching/matching_summary.csv" ]; then
    echo "[Step 8a] Running spectral library matching..."
    python3 "${SCRIPT_DIR}/scripts/matchms_match.py" \
        "${DATA_DIR}/spectral_library.msp" \
        "${OUTPUT_DIR}/camera/pseudospectra_info.csv" \
        "${OUTPUT_DIR}/camera/camera_peaklist.csv" \
        "${OUTPUT_DIR}/matching"
else
    echo "[Step 8a] Library matching already done, skipping."
fi

# Step 8b: Statistical summary is computed during report generation (parallel with 8a)

# --- CONVERGE2: matching and statistics merge into report ---

# ============================================================
# Step 9: Generate final report
# ============================================================
echo "[Step 9] Generating final report..."
python3 "${SCRIPT_DIR}/scripts/generate_report.py" \
    "${OUTPUT_DIR}" "${RESULTS_DIR}"

# ============================================================
# Validate report
# ============================================================
echo ""
echo "=== Validation ==="
if [ -f "${RESULTS_DIR}/report.csv" ]; then
    N_METRICS=$(tail -n +2 "${RESULTS_DIR}/report.csv" | wc -l)
    N_EMPTY=$(grep -cE ",,|,$" "${RESULTS_DIR}/report.csv" || true)
    echo "Report: ${N_METRICS} metrics, ${N_EMPTY} empty values"
    if [ "${N_EMPTY}" -gt 0 ]; then
        echo "WARNING: Report has empty values!"
    fi
    echo ""
    cat "${RESULTS_DIR}/report.csv"
else
    echo "ERROR: Report not generated!"
    exit 1
fi

echo ""
echo "=== Pipeline complete ==="
