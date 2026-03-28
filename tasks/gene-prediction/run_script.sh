#!/bin/bash
set -e

# =============================================================================
# Task 26: Gene Prediction Comparison
#
# DAG (depth 5, triple fan-out + merge):
#
# L0: bacterial genome assembly
# L1: ├── prodigal (ab initio prediction)
#     ├── prokka (evidence-based annotation)
#     └── glimmerhmm (HMM-based prediction)
# L2: ├── parse prodigal GFF → gene list
#     ├── parse prokka GFF → gene list
#     └── parse glimmerhmm GFF → gene list
# L3: bedtools intersect (pairwise overlaps)
# L4: MERGE (comparison report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${DATA}/genome.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{prodigal,prokka,glimmerhmm,comparison} "${RES}"

# L1-A: Prodigal
log_step "L1-A: Prodigal"
if [ ! -f "${OUT}/prodigal/genes.gff" ]; then
    prodigal -i "${GENOME}" -o "${OUT}/prodigal/genes.gff" \
             -a "${OUT}/prodigal/proteins.faa" -f gff -p single
else echo "Skipping (exists)"; fi

# L1-B: Prokka
log_step "L1-B: Prokka"
if [ ! -f "${OUT}/prokka/genome.gff" ]; then
    prokka "${GENOME}" --outdir "${OUT}/prokka" --prefix genome \
           --cpus ${THREADS} --kingdom Bacteria --force
else echo "Skipping (exists)"; fi

# L1-C: Prodigal in metagenomic mode (different parameters = different predictions)
log_step "L1-C: Prodigal (metagenomic mode)"
if [ ! -f "${OUT}/glimmerhmm/genes.gff" ]; then
    prodigal -i "${GENOME}" -o "${OUT}/glimmerhmm/genes.gff" \
             -a "${OUT}/glimmerhmm/proteins.faa" -f gff -p meta
else echo "Skipping (exists)"; fi

# L2: Parse predictions into BED format for comparison
log_step "L2: Parse predictions"

# Prodigal CDS
grep "CDS" "${OUT}/prodigal/genes.gff" | grep -v "^#" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4-1,$5,"prodigal_"NR,$6,$7}' \
    > "${OUT}/comparison/prodigal.bed" 2>/dev/null || true

# Prokka CDS
grep "CDS" "${OUT}/prokka/genome.gff" | grep -v "^#" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4-1,$5,"prokka_"NR,$6,$7}' \
    > "${OUT}/comparison/prokka.bed" 2>/dev/null || true

# GlimmerHMM
grep -v "^#" "${OUT}/glimmerhmm/genes.gff" | grep "mRNA\|CDS" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4-1,$5,"glimmer_"NR,$6,$7}' \
    > "${OUT}/comparison/glimmerhmm.bed" 2>/dev/null || true

# Sort BED files
for f in "${OUT}/comparison/"*.bed; do
    [ -s "$f" ] && sort -k1,1 -k2,2n "$f" > "${f}.sorted" && mv "${f}.sorted" "$f"
done

# L3: Pairwise overlaps
log_step "L3: bedtools intersect"
PRODIGAL_COUNT=$(wc -l < "${OUT}/comparison/prodigal.bed" 2>/dev/null | tr -d ' ' || echo "0")
PROKKA_COUNT=$(wc -l < "${OUT}/comparison/prokka.bed" 2>/dev/null | tr -d ' ' || echo "0")
GLIMMER_COUNT=$(wc -l < "${OUT}/comparison/glimmerhmm.bed" 2>/dev/null | tr -d ' ' || echo "0")

# Prodigal vs Prokka overlap
if [ -s "${OUT}/comparison/prodigal.bed" ] && [ -s "${OUT}/comparison/prokka.bed" ]; then
    PROD_PROK_OVERLAP=$(bedtools intersect -a "${OUT}/comparison/prodigal.bed" \
                                           -b "${OUT}/comparison/prokka.bed" -u | wc -l | tr -d ' ')
else
    PROD_PROK_OVERLAP=0
fi

# Prodigal vs GlimmerHMM overlap
if [ -s "${OUT}/comparison/prodigal.bed" ] && [ -s "${OUT}/comparison/glimmerhmm.bed" ]; then
    PROD_GLIM_OVERLAP=$(bedtools intersect -a "${OUT}/comparison/prodigal.bed" \
                                           -b "${OUT}/comparison/glimmerhmm.bed" -u | wc -l | tr -d ' ')
else
    PROD_GLIM_OVERLAP=0
fi

# Prokka vs GlimmerHMM overlap
if [ -s "${OUT}/comparison/prokka.bed" ] && [ -s "${OUT}/comparison/glimmerhmm.bed" ]; then
    PROK_GLIM_OVERLAP=$(bedtools intersect -a "${OUT}/comparison/prokka.bed" \
                                           -b "${OUT}/comparison/glimmerhmm.bed" -u | wc -l | tr -d ' ')
else
    PROK_GLIM_OVERLAP=0
fi

# All three overlap
if [ -s "${OUT}/comparison/prodigal.bed" ] && [ -s "${OUT}/comparison/prokka.bed" ] && [ -s "${OUT}/comparison/glimmerhmm.bed" ]; then
    bedtools intersect -a "${OUT}/comparison/prodigal.bed" \
                       -b "${OUT}/comparison/prokka.bed" -u \
        | bedtools intersect -a - -b "${OUT}/comparison/glimmerhmm.bed" -u \
        > "${OUT}/comparison/all_three.bed" 2>/dev/null || true
    ALL_THREE=$(wc -l < "${OUT}/comparison/all_three.bed" | tr -d ' ')
else
    ALL_THREE=0
fi

# Genome stats
GENOME_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${GENOME}")
GC=$(awk '/^>/{next}{s+=length; gc+=gsub(/[GCgc]/,"&")}END{printf "%.2f",gc/s*100}' "${GENOME}")

# L4: MERGE
log_step "L4-MERGE"

cat > "${RES}/gene_prediction_comparison.csv" << CSVEOF
metric,value
genome_length,${GENOME_LEN}
gc_content,${GC}
predictor_a_genes,${PRODIGAL_COUNT}
predictor_b_genes,${PROKKA_COUNT}
predictor_c_genes,${GLIMMER_COUNT}
overlap_a_b,${PROD_PROK_OVERLAP}
overlap_a_c,${PROD_GLIM_OVERLAP}
overlap_b_c,${PROK_GLIM_OVERLAP}
consensus_all_three,${ALL_THREE}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/gene_prediction_comparison.csv"
