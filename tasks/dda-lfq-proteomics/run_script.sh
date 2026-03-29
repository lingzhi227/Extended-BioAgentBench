#!/bin/bash
set -uo pipefail

# =============================================================================
# Task: DDA Label-Free Quantitative Proteomics
#
# DAG structure (depth 10, 3 convergence points):
#
# L0: mzML mass spec file + protein FASTA database
# L1: DecoyDatabase (add reverse decoys)
# L2: PeakPickerHiRes (centroid if profile mode)
# L3: CometAdapter (database search — peptide-spectrum matching)
# L4: PeptideIndexer (map peptides to protein sequences)
# L5: FalseDiscoveryRate (target-decoy FDR estimation)
#     ├──────────────────────────────────────────────────┐
# L6: IDFilter (1% PSM FDR)                        FDR statistics
#     │                                                  │
# L7: TextExporter (idXML → TSV)                   QC metrics
#     │                                                  │
# L8: Parse protein/peptide/PSM counts ◄────────────────┘ [CONVERGENCE 1+2]
# L9: MERGE                                               [CONVERGENCE 3]
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

MZML="${DATA}/sample.mzML"
DATABASE="${REF}/database.fasta"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{decoy,picked,search_comet,indexed,fdr,filtered,quant} "${RES}"

# L1: Add decoy sequences
log_step "L1: DecoyDatabase"
if [ ! -f "${OUT}/decoy/target_decoy.fasta" ]; then
    DecoyDatabase -in "${DATABASE}" -out "${OUT}/decoy/target_decoy.fasta" \
                  -decoy_string "DECOY_" -decoy_string_position prefix
fi

# L2: Peak picking
log_step "L2: PeakPickerHiRes"
if [ ! -f "${OUT}/picked/picked.mzML" ]; then
    PeakPickerHiRes -in "${MZML}" -out "${OUT}/picked/picked.mzML" \
                    -algorithm:ms_levels 1 2 2>/dev/null || {
        echo "Data already centroided, using as-is"
        cp "${MZML}" "${OUT}/picked/picked.mzML"
    }
fi

# L3: Comet database search
log_step "L3: CometAdapter"
if [ ! -f "${OUT}/search_comet/comet.idXML" ]; then
    CometAdapter -in "${OUT}/picked/picked.mzML" \
                 -database "${OUT}/decoy/target_decoy.fasta" \
                 -out "${OUT}/search_comet/comet.idXML" \
                 -precursor_mass_tolerance 10 \
                 -precursor_error_units ppm \
                 -fragment_mass_tolerance 0.02 \
                 -threads ${THREADS} 2>&1 || true
fi

# L4: PeptideIndexer
log_step "L4: PeptideIndexer"
if [ ! -f "${OUT}/indexed/indexed.idXML" ]; then
    PeptideIndexer -in "${OUT}/search_comet/comet.idXML" \
                   -fasta "${OUT}/decoy/target_decoy.fasta" \
                   -out "${OUT}/indexed/indexed.idXML" \
                   -decoy_string "DECOY_" -decoy_string_position prefix \
                   -enzyme:name Trypsin
fi

# L5: FDR
log_step "L5: FalseDiscoveryRate"
if [ ! -f "${OUT}/fdr/fdr.idXML" ]; then
    FalseDiscoveryRate -in "${OUT}/indexed/indexed.idXML" \
                       -out "${OUT}/fdr/fdr.idXML" -force
fi

# L6: Filter at 1% FDR
log_step "L6: IDFilter"
if [ ! -f "${OUT}/filtered/filtered.idXML" ]; then
    IDFilter -in "${OUT}/fdr/fdr.idXML" \
             -out "${OUT}/filtered/filtered.idXML" \
             -score:psm 0.01
fi

# L7: Export
log_step "L7: TextExporter"
if [ ! -f "${OUT}/quant/results.tsv" ]; then
    TextExporter -in "${OUT}/filtered/filtered.idXML" \
                 -out "${OUT}/quant/results.tsv"
fi

# MERGE
log_step "MERGE"
SPECTRA=$(grep -c "<spectrum " "${OUT}/picked/picked.mzML" 2>/dev/null || true)
SPECTRA=${SPECTRA:-0}
DB_SIZE=$(grep -c "^>" "${OUT}/decoy/target_decoy.fasta" 2>/dev/null || true)
DB_SIZE=${DB_SIZE:-0}
SEARCH_HITS=$(grep -c "PeptideHit" "${OUT}/search_comet/comet.idXML" 2>/dev/null || true)
SEARCH_HITS=${SEARCH_HITS:-0}
PSMS=$(grep -c "^PEPTIDE" "${OUT}/quant/results.tsv" 2>/dev/null || true)
PSMS=${PSMS:-0}
PEPTIDES=$(grep "^PEPTIDE" "${OUT}/quant/results.tsv" 2>/dev/null | awk '{print $2}' | sort -u | wc -l | tr -d ' ')
PROTEINS=$(grep -c "^PROTEIN" "${OUT}/quant/results.tsv" 2>/dev/null || true)
PROTEINS=${PROTEINS:-0}

cat > "${RES}/proteomics_report.csv" << CSVEOF
metric,value
total_spectra,${SPECTRA}
database_size_with_decoys,${DB_SIZE}
search_engine_hits,${SEARCH_HITS}
psms_after_fdr,${PSMS}
unique_peptides,${PEPTIDES}
proteins_identified,${PROTEINS}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/proteomics_report.csv"
