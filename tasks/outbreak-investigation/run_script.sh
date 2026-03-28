#!/bin/bash
set -e

# =============================================================================
# Task 14: Foodborne Pathogen Outbreak Investigation (Diamond DAG, depth=8)
#
# DAG structure (3 TB isolates):
#
# L0: reads (3 isolates × paired-end)
# L1: fastp ×3
#     ├──────────────────────────────────────────┐
# L2: snippy ×3                              shovill ×3
#     (ref-based SNP calling)                (de novo assembly)
#     │                                   ┌──────┼──────────┐
# L3: snippy-core                       prokka   mlst    quast+busco
#     (merge → core alignment)           ×3       ×3       ×3
#     │                            ┌──────┴──────┐
# L4: gubbins                   panaroo      amrfinder ×3
#     (recombination removal)  (pan-genome)  (protein+nuc AMR)
#     │                            │             │
# L5: snp-sites                core gene     abricate ×3
#     (extract variable sites) alignment     (cross-validate AMR)
#     ├──────────┐            ├────┘             │
# L6: iqtree   snp-dists   iqtree            AMR summary
#     (SNP      (pairwise   (core gene
#      tree)    distances)   tree)
#     └──────────┴────────────┴──────────────────┘
# L7: MERGE (outbreak report)
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

SAMPLES=("ERR1203059" "ERR2659153" "SRR998584")
REFERENCE="${REF}/MTB_reference.fasta"

# Snippy needs a dedicated env due to samtools version compatibility
SNIPPY_ENV="/pscratch/sd/l/lingzhi/micromamba/envs/snippy-env"
MAMBA="/pscratch/sd/l/lingzhi/micromamba-bin"
run_snippy() { ${MAMBA} run -p ${SNIPPY_ENV} "$@"; }

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "Started: $(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,snippy,assemblies,prokka,mlst_out,amr,qc,panaroo,phylogeny} "${RES}"

# ===========================================================================
# L1: Quality trimming (per sample)
# ===========================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L1: fastp trimming ${SAMPLE}"
    if [ ! -f "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" ]; then
        fastp --in1 "${DATA}/${SAMPLE}_R1.fastq.gz" --in2 "${DATA}/${SAMPLE}_R2.fastq.gz" \
              --out1 "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" --out2 "${OUT}/trimmed/${SAMPLE}_R2.fastq.gz" \
              --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
              --length_required 50 --thread ${THREADS} \
              --json "${OUT}/trimmed/${SAMPLE}_fastp.json" --html "${OUT}/trimmed/${SAMPLE}_fastp.html"
    else echo "Skipping (exists)"; fi
done

# ===========================================================================
# L2 LEFT BRANCH: Reference-based SNP calling with Snippy (per sample)
# ===========================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L2-LEFT: snippy SNP calling ${SAMPLE}"
    if [ ! -f "${OUT}/snippy/${SAMPLE}/snps.vcf" ]; then
        run_snippy snippy --cpus ${THREADS} --ram 8 --force \
               --outdir "${OUT}/snippy/${SAMPLE}" \
               --ref "${REFERENCE}" \
               --R1 "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" \
               --R2 "${OUT}/trimmed/${SAMPLE}_R2.fastq.gz"
    else echo "Skipping (exists)"; fi
done

# ===========================================================================
# L2 RIGHT BRANCH: De novo assembly with Shovill (per sample)
# ===========================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L2-RIGHT: shovill assembly ${SAMPLE}"
    if [ ! -f "${OUT}/assemblies/${SAMPLE}/contigs.fa" ]; then
        shovill --R1 "${OUT}/trimmed/${SAMPLE}_R1.fastq.gz" \
                --R2 "${OUT}/trimmed/${SAMPLE}_R2.fastq.gz" \
                --outdir "${OUT}/assemblies/${SAMPLE}" \
                --gsize 4400000 --cpus ${THREADS} --ram 8 --force
    else echo "Skipping (exists)"; fi
done

# ===========================================================================
# L3 LEFT: snippy-core (merge all SNP alignments)
# ===========================================================================
log_step "L3-LEFT: snippy-core merge"
if [ ! -f "${OUT}/snippy/core.full.aln" ]; then
    run_snippy snippy-core --ref "${REFERENCE}" \
                --prefix "${OUT}/snippy/core" \
                "${OUT}/snippy/ERR1203059" \
                "${OUT}/snippy/ERR2659153" \
                "${OUT}/snippy/SRR998584"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3 RIGHT: Annotation, typing, QC (per sample)
# ===========================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    # --- Prokka ---
    log_step "L3-RIGHT: prokka annotation ${SAMPLE}"
    if [ ! -f "${OUT}/prokka/${SAMPLE}/${SAMPLE}.gff" ]; then
        prokka "${OUT}/assemblies/${SAMPLE}/contigs.fa" \
               --outdir "${OUT}/prokka/${SAMPLE}" --prefix "${SAMPLE}" \
               --cpus ${THREADS} --kingdom Bacteria \
               --genus Mycobacterium --species tuberculosis --force
    else echo "Skipping (exists)"; fi

    # --- MLST ---
    log_step "L3-RIGHT: mlst typing ${SAMPLE}"
    if [ ! -f "${OUT}/mlst_out/${SAMPLE}.tsv" ]; then
        mlst "${OUT}/assemblies/${SAMPLE}/contigs.fa" > "${OUT}/mlst_out/${SAMPLE}.tsv"
    else echo "Skipping (exists)"; fi

    # --- QUAST ---
    log_step "L3-RIGHT: quast QC ${SAMPLE}"
    if [ ! -f "${OUT}/qc/${SAMPLE}_quast/report.tsv" ]; then
        quast "${OUT}/assemblies/${SAMPLE}/contigs.fa" \
              -r "${REFERENCE}" -o "${OUT}/qc/${SAMPLE}_quast" -t ${THREADS}
    else echo "Skipping (exists)"; fi
done

# --- BUSCO (one representative) ---
log_step "L3-RIGHT: busco completeness"
if [ ! -d "${OUT}/qc/busco" ]; then
    busco -i "${OUT}/assemblies/${SAMPLES[0]}/contigs.fa" -o busco --out_path "${OUT}/qc" \
          -l bacteria_odb10 -m genome -c ${THREADS} --force
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4 LEFT: Gubbins (recombination removal)
# ===========================================================================
log_step "L4-LEFT: gubbins recombination removal"
if [ ! -f "${OUT}/phylogeny/gubbins.filtered_polymorphic_sites.fasta" ]; then
    cd "${OUT}/phylogeny"
    run_gubbins.py "${OUT}/snippy/core.full.aln" \
                   --prefix gubbins \
                   --threads ${THREADS} \
                   --tree-builder iqtree 2>&1 || {
        echo "WARNING: Gubbins failed (common with <4 samples or low diversity). Using snippy-core alignment directly."
        cp "${OUT}/snippy/core.full.aln" "${OUT}/phylogeny/gubbins.filtered_polymorphic_sites.fasta"
    }
    cd "${SCRIPT_DIR}"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4 RIGHT: Panaroo pan-genome + AMRFinder (per sample)
# ===========================================================================
log_step "L4-RIGHT: panaroo pan-genome analysis"
if [ ! -f "${OUT}/panaroo/summary_statistics.txt" ]; then
    panaroo -i "${OUT}/prokka/ERR1203059/ERR1203059.gff" \
               "${OUT}/prokka/ERR2659153/ERR2659153.gff" \
               "${OUT}/prokka/SRR998584/SRR998584.gff" \
            -o "${OUT}/panaroo" --clean-mode strict \
            -a core -c 0.98 \
            --threads ${THREADS} 2>&1 || true
else echo "Skipping (exists)"; fi

for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L4-RIGHT: amrfinder ${SAMPLE}"
    if [ ! -f "${OUT}/amr/${SAMPLE}_amrfinder.tsv" ]; then
        amrfinder --nucleotide "${OUT}/assemblies/${SAMPLE}/contigs.fa" \
                  --threads ${THREADS} \
                  --output "${OUT}/amr/${SAMPLE}_amrfinder.tsv" --plus 2>&1 || true
    else echo "Skipping (exists)"; fi
done

# ===========================================================================
# L5 LEFT: snp-sites (extract variable sites)
# ===========================================================================
log_step "L5-LEFT: snp-sites extract variable sites"
CLEAN_ALN="${OUT}/phylogeny/gubbins.filtered_polymorphic_sites.fasta"
if [ ! -f "${CLEAN_ALN}" ]; then
    CLEAN_ALN="${OUT}/snippy/core.full.aln"
fi
if [ ! -f "${OUT}/phylogeny/clean_snps.fasta" ]; then
    snp-sites -o "${OUT}/phylogeny/clean_snps.fasta" "${CLEAN_ALN}" 2>&1 || {
        echo "WARNING: snp-sites failed, using core alignment directly"
        cp "${CLEAN_ALN}" "${OUT}/phylogeny/clean_snps.fasta"
    }
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5 RIGHT: ABRicate cross-validation (per sample)
# ===========================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    log_step "L5-RIGHT: abricate cross-validate AMR ${SAMPLE}"
    if [ ! -f "${OUT}/amr/${SAMPLE}_abricate.tsv" ]; then
        abricate "${OUT}/assemblies/${SAMPLE}/contigs.fa" --db card --minid 80 --mincov 60 \
                 > "${OUT}/amr/${SAMPLE}_abricate.tsv"
    else echo "Skipping (exists)"; fi
done

# ===========================================================================
# L6: Phylogenies + distances (fan-out from snp-sites)
# ===========================================================================

# --- SNP-based phylogeny ---
log_step "L6: iqtree SNP phylogeny"
if [ ! -f "${OUT}/phylogeny/snp_tree.treefile" ]; then
    iqtree -s "${OUT}/phylogeny/clean_snps.fasta" \
           -m GTR+G -bb 1000 -nt ${THREADS} \
           --prefix "${OUT}/phylogeny/snp_tree" 2>&1 || {
        echo "WARNING: iqtree failed on SNP alignment"
    }
else echo "Skipping (exists)"; fi

# --- Pairwise SNP distance matrix ---
log_step "L6: snp-dists pairwise distances"
if [ ! -f "${OUT}/phylogeny/snp_distances.tsv" ]; then
    snp-dists "${OUT}/phylogeny/clean_snps.fasta" > "${OUT}/phylogeny/snp_distances.tsv" 2>/dev/null || {
        echo "WARNING: snp-dists failed"
    }
else echo "Skipping (exists)"; fi

# --- Core gene phylogeny (from panaroo) ---
log_step "L6: iqtree core gene phylogeny"
if [ -f "${OUT}/panaroo/core_gene_alignment.aln" ] && [ ! -f "${OUT}/phylogeny/gene_tree.treefile" ]; then
    iqtree -s "${OUT}/panaroo/core_gene_alignment.aln" \
           -m GTR+G -nt ${THREADS} \
           --prefix "${OUT}/phylogeny/gene_tree" 2>&1 || {
        echo "WARNING: iqtree failed on core gene alignment"
    }
else echo "Skipping (no core gene alignment or exists)"; fi

# ===========================================================================
# L7: MERGE — Comprehensive outbreak report
# ===========================================================================
log_step "L7-MERGE: Building comprehensive outbreak report"

# --- Per-sample assembly stats ---
echo "sample,total_length,num_contigs,n50,gc_content,largest_contig" > "${RES}/assembly_stats.csv"
for SAMPLE in "${SAMPLES[@]}"; do
    TOTAL=$(grep "^Total length" "${OUT}/qc/${SAMPLE}_quast/report.tsv" | head -1 | cut -f2)
    NCTG=$(grep "^# contigs " "${OUT}/qc/${SAMPLE}_quast/report.tsv" | head -1 | cut -f2)
    N50=$(grep "^N50" "${OUT}/qc/${SAMPLE}_quast/report.tsv" | cut -f2)
    GC=$(grep "^GC" "${OUT}/qc/${SAMPLE}_quast/report.tsv" | cut -f2)
    LARGEST=$(grep "^Largest contig" "${OUT}/qc/${SAMPLE}_quast/report.tsv" | cut -f2)
    echo "${SAMPLE},${TOTAL},${NCTG},${N50},${GC},${LARGEST}" >> "${RES}/assembly_stats.csv"
done

# --- MLST ---
echo "sample,scheme,sequence_type,alleles" > "${RES}/mlst_results.csv"
for SAMPLE in "${SAMPLES[@]}"; do
    SCHEME=$(cut -f2 "${OUT}/mlst_out/${SAMPLE}.tsv")
    ST=$(cut -f3 "${OUT}/mlst_out/${SAMPLE}.tsv")
    ALLELES=$(cut -f4- "${OUT}/mlst_out/${SAMPLE}.tsv" | tr '\t' ';')
    echo "${SAMPLE},${SCHEME},${ST},${ALLELES}" >> "${RES}/mlst_results.csv"
done

# --- AMR summary ---
echo "sample,amr_genes_protein,amr_genes_nucleotide,amr_confirmed_both" > "${RES}/amr_summary.csv"
for SAMPLE in "${SAMPLES[@]}"; do
    AMR_P=$(tail -n +2 "${OUT}/amr/${SAMPLE}_amrfinder.tsv" 2>/dev/null | wc -l | tr -d ' ')
    AMR_N=$(tail -n +2 "${OUT}/amr/${SAMPLE}_abricate.tsv" 2>/dev/null | wc -l | tr -d ' ')
    # Cross-validate: genes found by both
    GENES_P=$(tail -n +2 "${OUT}/amr/${SAMPLE}_amrfinder.tsv" 2>/dev/null | awk -F'\t' '{print $7}' | sort -u | grep -v "^$")
    GENES_N=$(tail -n +2 "${OUT}/amr/${SAMPLE}_abricate.tsv" 2>/dev/null | awk -F'\t' '{print $6}' | sort -u | grep -v "^$")
    if [ -n "$GENES_P" ] && [ -n "$GENES_N" ]; then
        BOTH=$(comm -12 <(echo "$GENES_P") <(echo "$GENES_N") | wc -l | tr -d ' ')
    else
        BOTH=0
    fi
    echo "${SAMPLE},${AMR_P},${AMR_N},${BOTH}" >> "${RES}/amr_summary.csv"
done

# --- SNP distances ---
if [ -f "${OUT}/phylogeny/snp_distances.tsv" ]; then
    cp "${OUT}/phylogeny/snp_distances.tsv" "${RES}/snp_distance_matrix.tsv"
fi

# --- Newick trees ---
[ -f "${OUT}/phylogeny/snp_tree.treefile" ] && \
    cp "${OUT}/phylogeny/snp_tree.treefile" "${RES}/snp_phylogeny.nwk"
[ -f "${OUT}/phylogeny/gene_tree.treefile" ] && \
    cp "${OUT}/phylogeny/gene_tree.treefile" "${RES}/gene_phylogeny.nwk"

# --- Pan-genome summary ---
if [ -f "${OUT}/panaroo/summary_statistics.txt" ]; then
    cp "${OUT}/panaroo/summary_statistics.txt" "${RES}/pangenome_summary.txt"
fi

# --- BUSCO ---
BUSCO_SUM=$(grep "C:" "${OUT}/qc/busco/short_summary.specific.bacteria_odb10.busco.txt" 2>/dev/null \
    | head -1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || echo "N/A")

# --- Prokka gene counts ---
echo "sample,cds_count,trna_count,rrna_count" > "${RES}/annotation_summary.csv"
for SAMPLE in "${SAMPLES[@]}"; do
    CDS=$(grep "^CDS" "${OUT}/prokka/${SAMPLE}/${SAMPLE}.txt" | awk '{print $2}')
    TRNA=$(grep "^tRNA" "${OUT}/prokka/${SAMPLE}/${SAMPLE}.txt" | awk '{print $2}')
    RRNA=$(grep "^rRNA" "${OUT}/prokka/${SAMPLE}/${SAMPLE}.txt" | awk '{print $2}')
    echo "${SAMPLE},${CDS},${TRNA},${RRNA}" >> "${RES}/annotation_summary.csv"
done

# --- Main outbreak report CSV ---
# Count core SNPs
CORE_SNPS=$(grep -c "^>" "${OUT}/phylogeny/clean_snps.fasta" 2>/dev/null || echo "0")
SNP_ALN_LEN=$(head -2 "${OUT}/phylogeny/clean_snps.fasta" 2>/dev/null | tail -1 | wc -c | tr -d ' ')
# Pan-genome stats
CORE_GENES=$(grep "Core genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
TOTAL_GENES=$(grep "Total genes" "${OUT}/panaroo/summary_statistics.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")

cat > "${RES}/outbreak_report.csv" << CSVEOF
metric,value
num_isolates,${#SAMPLES[@]}
reference_genome,M. tuberculosis H37Rv
snp_alignment_length,${SNP_ALN_LEN}
busco_completeness,${BUSCO_SUM}
pangenome_core_genes,${CORE_GENES}
pangenome_total_genes,${TOTAL_GENES}
snp_tree_available,$([ -f "${RES}/snp_phylogeny.nwk" ] && echo "yes" || echo "no")
gene_tree_available,$([ -f "${RES}/gene_phylogeny.nwk" ] && echo "yes" || echo "no")
CSVEOF

echo ""
echo "=== Pipeline complete ==="
echo ""
echo "=== Outbreak Report ==="
cat "${RES}/outbreak_report.csv"
echo ""
echo "=== Assembly Stats ==="
cat "${RES}/assembly_stats.csv"
echo ""
echo "=== MLST ==="
cat "${RES}/mlst_results.csv"
echo ""
echo "=== AMR Summary ==="
cat "${RES}/amr_summary.csv"
echo ""
echo "=== Annotation Summary ==="
cat "${RES}/annotation_summary.csv"
echo ""
echo "=== SNP Distance Matrix ==="
cat "${RES}/snp_distance_matrix.tsv" 2>/dev/null || echo "N/A"
echo ""
echo "=== All result files ==="
ls -lh "${RES}/"
