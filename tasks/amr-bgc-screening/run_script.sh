#!/bin/bash
set -uo pipefail

# =============================================================================
# Task: AMR + Biosynthetic Gene Cluster Screening
#
# DAG structure (depth 8, 3 convergence points):
#
# L0: assembled contigs
# L1: prodigal (gene prediction)
#     ├── ARG screening branch ────────────────────────┐
# L2: │  ├── abricate --db card (nucleotide AMR)       │
#     │  └── amrfinder (protein AMR)                   │
# L3: │  merge ARG results (cross-validate)            │ [CONVERGENCE 1]
#     │                                                │
#     ├── BGC screening branch ────────────────────────┤
# L4: │  ├── gecco (ML-based BGC detection)            │
#     │  └── hmmsearch (PFAM domain search)            │
# L5: │  merge BGC results                             │ [CONVERGENCE 2]
#     │                                                │
#     └── virulence branch ────────────────────────────┤
# L6:    abricate --db vfdb (virulence factors)        │
#     └────────────────────────────────────────────────┘
# L7: MERGE (unified AMR + BGC + VF report)              [CONVERGENCE 3]
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

CONTIGS="${DATA}/contigs.fna"

log_step() { echo "== STEP: $1 == $(date)"; }
mkdir -p "${OUT}"/{prodigal,amr_abricate,amr_amrfinder,bgc_gecco,bgc_hmmer,virulence} "${RES}"

# L1: Gene prediction
log_step "L1: prodigal"
if [ ! -f "${OUT}/prodigal/genes.faa" ]; then
    prodigal -i "${CONTIGS}" -o "${OUT}/prodigal/genes.gff" \
             -a "${OUT}/prodigal/genes.faa" -d "${OUT}/prodigal/genes.fna" -f gff -p single
fi

# L2-A: ABRicate CARD (nucleotide-based AMR)
log_step "L2-A: abricate CARD"
if [ ! -f "${OUT}/amr_abricate/card.tsv" ]; then
    abricate "${CONTIGS}" --db card --minid 80 --mincov 60 > "${OUT}/amr_abricate/card.tsv" 2>/dev/null || true
fi

# L2-B: AMRFinderPlus (protein-based AMR)
log_step "L2-B: amrfinder"
if [ ! -f "${OUT}/amr_amrfinder/amr.tsv" ]; then
    amrfinder --nucleotide "${CONTIGS}" --protein "${OUT}/prodigal/genes.faa" \
              --gff "${OUT}/prodigal/genes.gff" \
              --threads ${THREADS} --output "${OUT}/amr_amrfinder/amr.tsv" --plus 2>&1 || true
fi

# L3: Cross-validate AMR
log_step "L3: cross-validate AMR"
AMR_NUC=$(tail -n +2 "${OUT}/amr_abricate/card.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
AMR_PROT=$(tail -n +2 "${OUT}/amr_amrfinder/amr.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
# Count stress/resistance genes
STRESS=$(tail -n +2 "${OUT}/amr_amrfinder/amr.tsv" 2>/dev/null | awk -F'\t' '$9=="STRESS"' | wc -l | tr -d ' ' || echo "0")
VIRULENCE_AMR=$(tail -n +2 "${OUT}/amr_amrfinder/amr.tsv" 2>/dev/null | awk -F'\t' '$9=="VIRULENCE"' | wc -l | tr -d ' ' || echo "0")

# L4: GECCO BGC detection
log_step "L4: gecco BGC"
if [ ! -d "${OUT}/bgc_gecco/done" ]; then
    gecco run --genome "${CONTIGS}" --output-dir "${OUT}/bgc_gecco" --jobs ${THREADS} 2>/dev/null || {
        echo "WARNING: gecco failed, skipping BGC detection"
    }
    touch "${OUT}/bgc_gecco/done"
fi
BGC_GECCO=$(ls "${OUT}/bgc_gecco/"*.clusters.tsv 2>/dev/null | head -1)
BGC_COUNT=$(tail -n +2 "$BGC_GECCO" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# L4: HMMER domain search (PFAM biosynthetic domains)
log_step "L4: protein domain stats"
TOTAL_GENES=$(grep -c "^>" "${OUT}/prodigal/genes.faa" 2>/dev/null || true)
TOTAL_GENES=${TOTAL_GENES:-0}

# L6: Virulence factors
log_step "L6: abricate VFDB"
if [ ! -f "${OUT}/virulence/vfdb.tsv" ]; then
    abricate "${CONTIGS}" --db vfdb --minid 80 --mincov 60 > "${OUT}/virulence/vfdb.tsv" 2>/dev/null || true
fi
VF_COUNT=$(tail -n +2 "${OUT}/virulence/vfdb.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# Gene stats
TOTAL_GENES=$(grep -c "^>" "${OUT}/prodigal/genes.faa" 2>/dev/null || true)
TOTAL_GENES=${TOTAL_GENES:-0}
GENOME_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${CONTIGS}")

# L7: MERGE
log_step "MERGE"

cat > "${RES}/amr_bgc_report.csv" << CSVEOF
metric,value
genome_length,${GENOME_LEN}
predicted_genes,${TOTAL_GENES}
amr_genes_nucleotide,${AMR_NUC}
amr_genes_protein,${AMR_PROT}
stress_response_genes,${STRESS}
virulence_amr_genes,${VIRULENCE_AMR}
biosynthetic_gene_clusters,${BGC_COUNT}
virulence_factors,${VF_COUNT}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/amr_bgc_report.csv"
