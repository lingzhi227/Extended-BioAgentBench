#!/bin/bash
set -e

# =============================================================================
# Task 19: Bacterial Pan-genome and Evolutionary Analysis
#
# DAG (depth 6, fan-out from annotation):
#
# L0: 5 E. coli genome assemblies
# L1: Prokka annotation (×5 parallel)
# L2: ├── Panaroo (pan-genome from 5 GFFs)
#     │     ├── L3: core gene alignment
#     │     │        └── L4: snp-sites → iqtree (core phylogeny)
#     │     └── L3: gene presence/absence matrix
#     │              └── L4: accessory genome stats
#     └── mlst (×5)
#     └── abricate (×5)
# L5: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

SAMPLES=("ecoli_K12" "ecoli_O157H7" "ecoli_CFT073" "ecoli_UTI89" "ecoli_APEC")

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{prokka,mlst_out,amr,panaroo,phylogeny} "${RES}"

# ===========================================================================
# L0: Decompress genomes if needed
# ===========================================================================
log_step "L0: Prepare genomes"
for SAMPLE in "${SAMPLES[@]}"; do
    if [ ! -f "${DATA}/${SAMPLE}.fna" ]; then
        gunzip -k "${DATA}/${SAMPLE}.fna.gz" 2>/dev/null || true
    fi
done

# ===========================================================================
# L1: Prokka annotation (×5)
# ===========================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L1: Prokka ${SAMPLE}"
    if [ ! -f "${OUT}/prokka/${SAMPLE}/${SAMPLE}.gff" ]; then
        prokka "${DATA}/${SAMPLE}.fna" --outdir "${OUT}/prokka/${SAMPLE}" --prefix "${SAMPLE}" \
               --cpus ${THREADS} --kingdom Bacteria --genus Escherichia --species coli --force
    else echo "Skipping (exists)"; fi
done

# ===========================================================================
# L1 parallel: MLST + AMR (×5)
# ===========================================================================
echo "sample,scheme,sequence_type" > "${OUT}/mlst_out/mlst_all.csv"
for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L1: mlst + abricate ${SAMPLE}"
    mlst "${DATA}/${SAMPLE}.fna" > "${OUT}/mlst_out/${SAMPLE}.tsv" 2>/dev/null || true
    SCHEME=$(cut -f2 "${OUT}/mlst_out/${SAMPLE}.tsv" 2>/dev/null)
    ST=$(cut -f3 "${OUT}/mlst_out/${SAMPLE}.tsv" 2>/dev/null)
    echo "${SAMPLE},${SCHEME},${ST}" >> "${OUT}/mlst_out/mlst_all.csv"

    if [ ! -f "${OUT}/amr/${SAMPLE}.tsv" ]; then
        abricate "${DATA}/${SAMPLE}.fna" --db card --minid 80 --mincov 60 > "${OUT}/amr/${SAMPLE}.tsv" 2>/dev/null || true
    fi
done

# ===========================================================================
# L2: Panaroo pan-genome
# ===========================================================================
log_step "L2: Panaroo"
if [ ! -f "${OUT}/panaroo/summary_statistics.txt" ]; then
    GFF_FILES=""
    for SAMPLE in "${SAMPLES[@]}"; do
        GFF_FILES="${GFF_FILES} ${OUT}/prokka/${SAMPLE}/${SAMPLE}.gff"
    done
    panaroo -i ${GFF_FILES} -o "${OUT}/panaroo" --clean-mode strict \
            -a core -c 0.98 --threads ${THREADS} 2>&1 || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3-L4: Core gene phylogeny
# ===========================================================================
log_step "L3: snp-sites on core alignment"
if [ -f "${OUT}/panaroo/core_gene_alignment.aln" ] && [ ! -f "${OUT}/phylogeny/core_snps.fasta" ]; then
    snp-sites -o "${OUT}/phylogeny/core_snps.fasta" "${OUT}/panaroo/core_gene_alignment.aln" 2>&1 || \
        cp "${OUT}/panaroo/core_gene_alignment.aln" "${OUT}/phylogeny/core_snps.fasta"
fi

log_step "L4: iqtree phylogeny"
if [ -f "${OUT}/phylogeny/core_snps.fasta" ] && [ ! -f "${OUT}/phylogeny/core_tree.treefile" ]; then
    iqtree -s "${OUT}/phylogeny/core_snps.fasta" -m GTR+G -bb 1000 \
           -nt ${THREADS} --prefix "${OUT}/phylogeny/core_tree" 2>&1 || true
fi

# ===========================================================================
# L5: MERGE
# ===========================================================================
log_step "L5-MERGE"

CORE_GENES=$(grep "Core genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
TOTAL_GENES=$(grep "Total genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
SHELL_GENES=$(grep "Shell genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
CLOUD_GENES=$(grep "Cloud genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")

cat > "${RES}/pangenome_report.csv" << CSVEOF
metric,value
num_genomes,${#SAMPLES[@]}
core_genes,${CORE_GENES}
shell_genes,${SHELL_GENES}
cloud_genes,${CLOUD_GENES}
total_genes,${TOTAL_GENES}
tree_available,$([ -f "${OUT}/phylogeny/core_tree.treefile" ] && echo "yes" || echo "no")
CSVEOF

cp "${OUT}/mlst_out/mlst_all.csv" "${RES}/" 2>/dev/null || true
cp "${OUT}/phylogeny/core_tree.treefile" "${RES}/core_phylogeny.nwk" 2>/dev/null || true
cp "${OUT}/panaroo/gene_presence_absence.csv" "${RES}/" 2>/dev/null || true

# Per-genome AMR summary
echo "sample,amr_genes" > "${RES}/amr_summary.csv"
for SAMPLE in "${SAMPLES[@]}"; do
    COUNT=$(tail -n +2 "${OUT}/amr/${SAMPLE}.tsv" 2>/dev/null | wc -l | tr -d ' ')
    echo "${SAMPLE},${COUNT}" >> "${RES}/amr_summary.csv"
done

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/pangenome_report.csv"
echo ""
cat "${RES}/mlst_all.csv"
echo ""
ls -lh "${RES}/"
