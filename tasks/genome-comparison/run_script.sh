#!/bin/bash
set -e

# =============================================================================
# Task 22: Pairwise Bacterial Genome Comparison
#
# DAG (depth 5):
#
# L0: genome_A.fna, genome_B.fna
# L1: ├── prokka A          prokka B
#     └── nucmer (whole-genome alignment A vs B)
# L2: ├── dnadiff (alignment summary)
#     │     └── show-snps (SNP extraction)
#     └── panaroo (ortholog clustering from prokka GFFs)
# L3: ├── SNP annotation (bedtools intersect SNPs with genes)
#     └── unique gene extraction per genome
# L4: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME_A="${DATA}/genome_A.fna"
GENOME_B="${DATA}/genome_B.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{prokka_A,prokka_B,nucmer,panaroo,snp_annotation} "${RES}"

# ===========================================================================
# L1-A: Annotate genome A
# ===========================================================================
log_step "L1-A: Prokka genome A"
if [ ! -f "${OUT}/prokka_A/genome_A.gff" ]; then
    prokka "${GENOME_A}" --outdir "${OUT}/prokka_A" --prefix genome_A \
           --cpus ${THREADS} --kingdom Bacteria --force
else echo "Skipping (exists)"; fi

# ===========================================================================
# L1-B: Annotate genome B
# ===========================================================================
log_step "L1-B: Prokka genome B"
if [ ! -f "${OUT}/prokka_B/genome_B.gff" ]; then
    prokka "${GENOME_B}" --outdir "${OUT}/prokka_B" --prefix genome_B \
           --cpus ${THREADS} --kingdom Bacteria --force
else echo "Skipping (exists)"; fi

# ===========================================================================
# L1-C: Whole-genome alignment with nucmer
# ===========================================================================
log_step "L1-C: nucmer alignment"
if [ ! -f "${OUT}/nucmer/alignment.delta" ]; then
    nucmer --prefix="${OUT}/nucmer/alignment" "${GENOME_A}" "${GENOME_B}"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2-A: dnadiff summary statistics
# ===========================================================================
log_step "L2-A: dnadiff"
if [ ! -f "${OUT}/nucmer/dnadiff.report" ]; then
    dnadiff -d "${OUT}/nucmer/alignment.delta" -p "${OUT}/nucmer/dnadiff"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2-B: Extract SNPs
# ===========================================================================
log_step "L2-B: show-snps"
if [ ! -f "${OUT}/nucmer/snps.tsv" ]; then
    show-snps -Clr "${OUT}/nucmer/alignment.delta" > "${OUT}/nucmer/snps.tsv"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2-C: Panaroo ortholog comparison
# ===========================================================================
log_step "L2-C: Panaroo ortholog comparison"
if [ ! -f "${OUT}/panaroo/summary_statistics.txt" ]; then
    panaroo -i "${OUT}/prokka_A/genome_A.gff" "${OUT}/prokka_B/genome_B.gff" \
            -o "${OUT}/panaroo" --clean-mode strict -c 0.98 \
            --threads ${THREADS} 2>&1 || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3: Parse results
# ===========================================================================
log_step "L3: Parse results"

# dnadiff stats
ALIGNED_BASES=$(grep "AlignedBases" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}' | sed 's/(.*//')
AVG_IDENTITY=$(grep "AvgIdentity" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}')
TOTAL_SNPS=$(grep "TotalSNPs" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}')
TOTAL_INDELS=$(grep "TotalIndels" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}')
BREAKPOINTS=$(grep "Breakpoints" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}')
RELOCATIONS=$(grep "Relocations" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}')
INVERSIONS=$(grep "Inversions" "${OUT}/nucmer/dnadiff.report" | head -1 | awk '{print $2}')

# Genome sizes
LEN_A=$(awk '/^>/{if(l)s+=l; l=0; next}{l+=length}END{s+=l; print s}' "${GENOME_A}")
LEN_B=$(awk '/^>/{if(l)s+=l; l=0; next}{l+=length}END{s+=l; print s}' "${GENOME_B}")

# Gene counts from prokka
CDS_A=$(grep "^CDS" "${OUT}/prokka_A/genome_A.txt" | awk '{print $2}')
CDS_B=$(grep "^CDS" "${OUT}/prokka_B/genome_B.txt" | awk '{print $2}')

# Panaroo: shared vs unique genes
CORE=$(grep "Core genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "0")
TOTAL_PAN=$(grep "Total genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "0")
UNIQUE=$((TOTAL_PAN - CORE))

# ===========================================================================
# L4: MERGE
# ===========================================================================
log_step "L4-MERGE"

cat > "${RES}/genome_comparison.csv" << CSVEOF
metric,value
genome_a_length,${LEN_A}
genome_b_length,${LEN_B}
aligned_bases,${ALIGNED_BASES}
average_identity,${AVG_IDENTITY}
total_snps,${TOTAL_SNPS}
total_indels,${TOTAL_INDELS}
breakpoints,${BREAKPOINTS}
relocations,${RELOCATIONS}
inversions,${INVERSIONS}
cds_genome_a,${CDS_A}
cds_genome_b,${CDS_B}
shared_orthologs,${CORE}
unique_genes,${UNIQUE}
total_pangenome,${TOTAL_PAN}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/genome_comparison.csv"
echo ""
ls -lh "${RES}/"
