#!/usr/bin/env bash
set -euo pipefail
#
# DDA-LFQ Proteomics Pipeline
#
# DAG Structure (depth=8, convergence=2):
#
#   mzML_1 → PeakPicker(1) → SearchEngine(2) → PeptideIndexer(3) ─┐
#   mzML_2 → PeakPicker(1) → SearchEngine(2) → PeptideIndexer(3) ─┤
#   ...                                                              │
#   mzML_N → PeakPicker(1) → SearchEngine(2) → PeptideIndexer(3) ─┤ CONVERGE1
#                                                                    ↓
#     PSMFeatureExtractor(4) → PercolatorAdapter(5) → IDFilter(6)
#       ├→ ProteinQuantifier(7a) ──┐
#       └→ TextExporter(7b) ───────┤ CONVERGE2
#                                   ↓
#                             report(8)
#
# Inputs:  data/*.mzML, reference/database.fasta, reference/design.tsv
# Outputs: results/report.csv

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
REF_DIR="${SCRIPT_DIR}/reference"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"
FASTA="${REF_DIR}/database.fasta"

mkdir -p "${OUTPUT_DIR}"/{picked,search,indexed,features,percolator,filtered,quant} "${RESULTS_DIR}"

echo "=== DDA-LFQ Proteomics Pipeline ==="
echo "Threads: ${THREADS}"

MZML_FILES=$(ls "${DATA_DIR}"/*.mzML)
N_RUNS=$(echo "$MZML_FILES" | wc -l)
echo "Found ${N_RUNS} mzML files"

for mzml in ${MZML_FILES}; do
    BASE=$(basename "$mzml" .mzML)

    # Step 1: Peak picking
    if [ ! -f "${OUTPUT_DIR}/picked/${BASE}.mzML" ]; then
        echo "[Step 1] PeakPicker: ${BASE}..."
        PeakPickerHiRes \
            -in "$mzml" \
            -out "${OUTPUT_DIR}/picked/${BASE}.mzML" \
            -threads ${THREADS} 2>&1 | tail -2
    fi

    # Step 2: Database search (Comet)
    if [ ! -f "${OUTPUT_DIR}/search/${BASE}.idXML" ]; then
        echo "[Step 2] Database search: ${BASE}..."
        CometAdapter \
            -in "${OUTPUT_DIR}/picked/${BASE}.mzML" \
            -out "${OUTPUT_DIR}/search/${BASE}.idXML" \
            -database "${FASTA}" \
            -comet_executable comet \
            -precursor_mass_tolerance 10 \
            -fragment_mass_tolerance 0.02 \
            -threads ${THREADS} 2>&1 | tail -2
    fi

    # Step 3: PeptideIndexer (add protein context)
    if [ ! -f "${OUTPUT_DIR}/indexed/${BASE}.idXML" ]; then
        echo "[Step 3] PeptideIndexer: ${BASE}..."
        PeptideIndexer \
            -in "${OUTPUT_DIR}/search/${BASE}.idXML" \
            -out "${OUTPUT_DIR}/indexed/${BASE}.idXML" \
            -fasta "${FASTA}" \
            -decoy_string _rev \
            -decoy_string_position suffix \
            -enzyme:specificity none 2>&1 | tail -2 || true
    fi

    # Step 4: PSMFeatureExtractor
    if [ ! -f "${OUTPUT_DIR}/features/${BASE}.idXML" ]; then
        echo "[Step 4] Feature extraction: ${BASE}..."
        PSMFeatureExtractor \
            -in "${OUTPUT_DIR}/indexed/${BASE}.idXML" \
            -out "${OUTPUT_DIR}/features/${BASE}.idXML" 2>&1 | tail -2
    fi
done

# Step 5: Merge all with merge_proteins_add_psms (CONVERGE1)
MERGED="${OUTPUT_DIR}/percolator/merged_all.idXML"
if [ ! -f "${MERGED}" ]; then
    echo "[Step 5a] Merging all search results..."
    MERGE_IN=""
    for mzml in ${MZML_FILES}; do
        BASE=$(basename "$mzml" .mzML)
        MERGE_IN="${MERGE_IN} ${OUTPUT_DIR}/features/${BASE}.idXML"
    done
    IDMerger -in ${MERGE_IN} -out "${MERGED}" \
        -merge_proteins_add_PSMs 2>&1 | tail -2
fi

# Step 5b: Percolator FDR scoring
SCORED="${OUTPUT_DIR}/percolator/scored.idXML"
if [ ! -f "${SCORED}" ]; then
    echo "[Step 5b] Percolator FDR scoring..."
    PercolatorAdapter \
        -in "${MERGED}" \
        -out "${SCORED}" \
        -percolator_executable percolator \
        -decoy_pattern _rev \
        -enzyme no_enzyme \
        -threads ${THREADS} 2>&1 | tail -5
fi

# Step 6: IDFilter (apply FDR threshold)
FILTERED="${OUTPUT_DIR}/filtered/filtered.idXML"
if [ ! -f "${FILTERED}" ]; then
    echo "[Step 6] IDFilter..."
    IDFilter \
        -in "${SCORED}" \
        -out "${FILTERED}" \
        -score:psm 0.05 2>&1 | tail -2
fi

# Step 7: Quantification and report (CONVERGE2)
echo "[Step 7-8] Generating report..."
python3 "${SCRIPT_DIR}/scripts/generate_report.py" \
    "${OUTPUT_DIR}" "${RESULTS_DIR}" "${N_RUNS}"

echo ""
echo "=== Final Report ==="
cat "${RESULTS_DIR}/report.csv"
echo ""
echo "=== Pipeline complete ==="
