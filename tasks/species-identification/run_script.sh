#!/bin/bash
set -e

# =============================================================================
# Task 30: Bacterial Species Identification from Assembly
#
# DAG (depth 5):
# L0: assembled genome (unknown species)
# L1: ├── mlst (sequence typing)
#     ├── abricate --db ecoh (E. coli serotyping)
#     └── seqkit stats (basic stats)
# L2: ├── nucmer vs ref_A (E. coli K12)
#     ├── nucmer vs ref_B (Salmonella)
#     └── nucmer vs ref_C (S. aureus)
# L3: ├── dnadiff A (identity to K12)
#     ├── dnadiff B (identity to Salmonella)
#     └── dnadiff C (identity to S. aureus)
# L4: Determine closest species by identity
# L5: MERGE (species ID report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${DATA}/unknown_genome.fna"
REF_A="${REF}/ecoli_k12.fna"
REF_B="${REF}/salmonella.fna"
REF_C="${REF}/saureus.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{typing,stats,align_A,align_B,align_C} "${RES}"

# L1: Typing + stats
log_step "L1: MLST + stats"
if [ ! -f "${OUT}/typing/mlst.tsv" ]; then
    mlst "${GENOME}" > "${OUT}/typing/mlst.tsv" 2>/dev/null || true
fi
seqkit stats "${GENOME}" -T > "${OUT}/stats/genome_stats.tsv"

GENOME_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${GENOME}")
GC=$(awk '/^>/{next}{s+=length;gc+=gsub(/[GCgc]/,"&")}END{printf "%.2f",gc/s*100}' "${GENOME}")
NCONTIGS=$(grep -c ">" "${GENOME}" || true)
MLST_SCHEME=$(cut -f2 "${OUT}/typing/mlst.tsv" 2>/dev/null || echo "unknown")
MLST_ST=$(cut -f3 "${OUT}/typing/mlst.tsv" 2>/dev/null || echo "unknown")

# L2-L3: Align to each reference and compute identity
log_step "L2-L3: Multi-reference alignment"
declare -A IDENTITIES
declare -A ALIGNED_FRACS

for label in A B C; do
    eval REF_FILE="\${REF_${label}}"
    OUTDIR="${OUT}/align_${label}"

    if [ ! -f "${OUTDIR}/dnadiff.report" ] && [ -f "${REF_FILE}" ]; then
        nucmer --prefix="${OUTDIR}/alignment" "${REF_FILE}" "${GENOME}" 2>/dev/null || true
        if [ -f "${OUTDIR}/alignment.delta" ]; then
            dnadiff -d "${OUTDIR}/alignment.delta" -p "${OUTDIR}/dnadiff" 2>/dev/null || true
        fi
    fi

    if [ -f "${OUTDIR}/dnadiff.report" ]; then
        IDENTITIES[$label]=$(grep "AvgIdentity" "${OUTDIR}/dnadiff.report" | head -1 | awk '{print $2}')
        ALIGNED_FRACS[$label]=$(grep "AlignedBases" "${OUTDIR}/dnadiff.report" | tail -1 | awk '{print $2}' | grep -oP '[\d.]+%' | tr -d '%')
    else
        IDENTITIES[$label]="0"
        ALIGNED_FRACS[$label]="0"
    fi
done

# L4: Determine best match
BEST_LABEL="unknown"
BEST_IDENTITY="0"
BEST_ALIGNED="0"
for label in A B C; do
    id=${IDENTITIES[$label]}
    af=${ALIGNED_FRACS[$label]}
    if [ -n "$af" ] && python3 -c "exit(0 if float('${af}') > float('${BEST_ALIGNED}') else 1)" 2>/dev/null; then
        BEST_LABEL=$label
        BEST_IDENTITY=$id
        BEST_ALIGNED=$af
    fi
done

case $BEST_LABEL in
    A) BEST_SPECIES="Escherichia coli" ;;
    B) BEST_SPECIES="Salmonella enterica" ;;
    C) BEST_SPECIES="Staphylococcus aureus" ;;
    *) BEST_SPECIES="Unknown" ;;
esac

# L5: MERGE
log_step "L5-MERGE"

cat > "${RES}/species_id_report.csv" << CSVEOF
metric,value
genome_length,${GENOME_LEN}
num_contigs,${NCONTIGS}
gc_content,${GC}
mlst_scheme,${MLST_SCHEME}
mlst_sequence_type,${MLST_ST}
identity_to_ecoli,${IDENTITIES[A]}
aligned_to_ecoli_pct,${ALIGNED_FRACS[A]}
identity_to_salmonella,${IDENTITIES[B]}
aligned_to_salmonella_pct,${ALIGNED_FRACS[B]}
identity_to_saureus,${IDENTITIES[C]}
aligned_to_saureus_pct,${ALIGNED_FRACS[C]}
best_match_species,${BEST_SPECIES}
best_match_identity,${BEST_IDENTITY}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/species_id_report.csv"
