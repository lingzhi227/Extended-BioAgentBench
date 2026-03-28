#!/bin/bash
set -e

# =============================================================================
# Task 13: Mobile Genetic Element Characterization (Diamond DAG, depth=5)
#
# DAG structure:
#
#   L0: reads
#   L1: fastp (trim)
#   L2: shovill (assemble)
#       │
#   L3: ├── prokka (annotate) ──────────────────────────────────────┐
#       │      ├── L4a: isescan (IS elements)                       │
#       │      └── L4b: integron_finder (integrons)                 │
#       │                                                           │
#       ├── amrfinder (protein AMR) ───┐                            │
#       │                              ├─ L4c: cross_validate_amr   │
#       ├── abricate CARD (nuc AMR) ───┘                            │
#       │                                                           │
#       ├── mob_recon (plasmid recon) ─── L4d: plasmidfinder.py ────┤
#       │                                 (replicon typing)         │
#       ├── abricate VFDB (virulence) ─────────────────────────────┤
#       │                                                           │
#       └── quast + busco (QC) ────────────────────────────────────┤
#                                                                   │
#   L5: MERGE ──────────────────────────────────────────────────────┘
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "Started: $(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,assembly,prokka,isescan,integrons,amr,plasmids,virulence,qc} "${RES}"

# ===========================================================================
# L0-L1: Preprocessing
# ===========================================================================
log_step "L1: Quality trimming with fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq" --in2 "${DATA}/reads_R2.fastq" \
          --out1 "${OUT}/trimmed/R1.fastq" --out2 "${OUT}/trimmed/R2.fastq" \
          --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
          --length_required 30 --thread ${THREADS} \
          --json "${OUT}/trimmed/fastp.json" --html "${OUT}/trimmed/fastp.html"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2: Assembly
# ===========================================================================
log_step "L2: Genome assembly with shovill"
if [ ! -f "${OUT}/assembly/contigs.fa" ]; then
    shovill --R1 "${OUT}/trimmed/R1.fastq" --R2 "${OUT}/trimmed/R2.fastq" \
            --outdir "${OUT}/assembly" --gsize 2914567 --cpus ${THREADS} --ram 8 --force
else echo "Skipping (exists)"; fi
CONTIGS="${OUT}/assembly/contigs.fa"

# ===========================================================================
# L3: DIAMOND BRANCHES (all depend only on contigs.fa)
# ===========================================================================

# --- L3-A: Gene annotation with Prokka ---
log_step "L3-A: Gene annotation with Prokka"
if [ ! -f "${OUT}/prokka/MRSA.gff" ]; then
    prokka "${CONTIGS}" --outdir "${OUT}/prokka" --prefix MRSA \
           --cpus ${THREADS} --kingdom Bacteria --genus Staphylococcus --species aureus --force
else echo "Skipping (exists)"; fi

# --- L3-B: AMR detection with AMRFinderPlus (protein-homology method) ---
log_step "L3-B: AMR detection with AMRFinderPlus"
if [ ! -f "${OUT}/amr/amrfinder_results.tsv" ]; then
    amrfinder --nucleotide "${CONTIGS}" --organism Staphylococcus_aureus \
              --threads ${THREADS} --output "${OUT}/amr/amrfinder_results.tsv" --plus 2>&1 || true
else echo "Skipping (exists)"; fi

# --- L3-C: AMR detection with ABRicate CARD (nucleotide method) ---
log_step "L3-C: AMR detection with ABRicate (CARD)"
if [ ! -f "${OUT}/amr/abricate_card.tsv" ]; then
    abricate "${CONTIGS}" --db card --minid 80 --mincov 60 > "${OUT}/amr/abricate_card.tsv"
else echo "Skipping (exists)"; fi

# --- L3-D: Plasmid reconstruction with mob_recon ---
log_step "L3-D: Plasmid reconstruction with mob_recon"
if [ ! -f "${OUT}/plasmids/contig_report.txt" ]; then
    mob_recon --infile "${CONTIGS}" --outdir "${OUT}/plasmids" \
              --num_threads ${THREADS} --force 2>&1 || true
else echo "Skipping (exists)"; fi

# --- L3-E: Virulence factor detection with ABRicate VFDB ---
log_step "L3-E: Virulence detection with ABRicate (VFDB)"
if [ ! -f "${OUT}/virulence/abricate_vfdb.tsv" ]; then
    abricate "${CONTIGS}" --db vfdb --minid 80 --mincov 60 > "${OUT}/virulence/abricate_vfdb.tsv"
else echo "Skipping (exists)"; fi

# --- L3-F: Assembly QC with QUAST ---
log_step "L3-F: Assembly QC with QUAST"
if [ ! -f "${OUT}/qc/quast/report.tsv" ]; then
    quast "${CONTIGS}" -o "${OUT}/qc/quast" -t ${THREADS}
else echo "Skipping (exists)"; fi

# --- L3-G: Assembly completeness with BUSCO ---
log_step "L3-G: Assembly completeness with BUSCO"
if [ ! -d "${OUT}/qc/busco" ]; then
    busco -i "${CONTIGS}" -o busco --out_path "${OUT}/qc" \
          -l bacteria_odb10 -m genome -c ${THREADS} --force
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4: SUB-BRANCHES (depend on L3 outputs → deeper diamond)
# ===========================================================================

# --- L4-A: IS element detection with ISEScan (depends on Prokka proteome/contigs) ---
log_step "L4-A: IS element detection with ISEScan"
if [ ! -d "${OUT}/isescan/hmm" ]; then
    cp "${CONTIGS}" "${OUT}/isescan/MRSA.fa"
    isescan.py --seqfile "${OUT}/isescan/MRSA.fa" --output "${OUT}/isescan" --nthread ${THREADS} || true
else echo "Skipping (exists)"; fi

# --- L4-B: Integron detection with integron_finder (depends on annotation context) ---
log_step "L4-B: Integron detection with integron_finder"
if [ ! -d "${OUT}/integrons/Results_Integron_Finder_MRSA" ]; then
    cp "${CONTIGS}" "${OUT}/integrons/MRSA.fasta"
    integron_finder --local-max --cpu ${THREADS} \
                    --outdir "${OUT}/integrons" \
                    "${OUT}/integrons/MRSA.fasta" 2>&1 || true
else echo "Skipping (exists)"; fi

# --- L4-C: Cross-validate AMR (merge AMRFinderPlus + ABRicate CARD results) ---
log_step "L4-C: Cross-validate AMR results"
echo "gene,method,identity" > "${OUT}/amr/cross_validated.csv"
# Genes found by AMRFinderPlus
if [ -s "${OUT}/amr/amrfinder_results.tsv" ]; then
    tail -n +2 "${OUT}/amr/amrfinder_results.tsv" | \
        awk -F'\t' 'BEGIN{OFS=","} {print $7,"amrfinder",$3}' >> "${OUT}/amr/cross_validated.csv"
fi
# Genes found by ABRicate CARD
if [ -s "${OUT}/amr/abricate_card.tsv" ]; then
    tail -n +2 "${OUT}/amr/abricate_card.tsv" | \
        awk -F'\t' 'BEGIN{OFS=","} {print $6,"abricate_card",$10}' >> "${OUT}/amr/cross_validated.csv"
fi
AMR_BOTH=$(awk -F',' 'NR>1{genes[$1]++} END{c=0; for(g in genes) if(genes[g]>=2) c++; print c}' \
    "${OUT}/amr/cross_validated.csv")
echo "AMR genes confirmed by both methods: ${AMR_BOTH}"

# --- L4-D: Replicon typing with plasmidfinder (depends on mob_recon output) ---
log_step "L4-D: Replicon typing with plasmidfinder"
if [ ! -f "${OUT}/plasmids/plasmidfinder_results.tsv" ]; then
    # Run plasmidfinder on each plasmid contig identified by mob_recon
    PLASMID_FA="${OUT}/plasmids/plasmid_contigs.fasta"
    if [ -f "${OUT}/plasmids/plasmid_AA109.fasta" ] || ls "${OUT}/plasmids/plasmid_"*.fasta 1>/dev/null 2>&1; then
        cat "${OUT}/plasmids/plasmid_"*.fasta > "${PLASMID_FA}" 2>/dev/null || true
    fi
    if [ -s "${PLASMID_FA}" ]; then
        plasmidfinder.py -i "${PLASMID_FA}" -o "${OUT}/plasmids/pf_out" -x 2>&1 || true
        cp "${OUT}/plasmids/pf_out/results_tab.tsv" "${OUT}/plasmids/plasmidfinder_results.tsv" 2>/dev/null || true
    else
        # No plasmid contigs found by mob_recon — run on full assembly
        plasmidfinder.py -i "${CONTIGS}" -o "${OUT}/plasmids/pf_out" -x 2>&1 || true
        cp "${OUT}/plasmids/pf_out/results_tab.tsv" "${OUT}/plasmids/plasmidfinder_results.tsv" 2>/dev/null || true
    fi
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5: MERGE — Combine all branch results
# ===========================================================================
log_step "L5-MERGE: Building comprehensive results CSV"

# Extract QUAST
TOTAL_LEN=$(grep "^Total length" "${OUT}/qc/quast/report.tsv" | head -1 | cut -f2)
NUM_CONTIGS=$(grep "^# contigs " "${OUT}/qc/quast/report.tsv" | head -1 | cut -f2)
N50=$(grep "^N50" "${OUT}/qc/quast/report.tsv" | cut -f2)
GC=$(grep "^GC" "${OUT}/qc/quast/report.tsv" | cut -f2)
LARGEST=$(grep "^Largest contig" "${OUT}/qc/quast/report.tsv" | cut -f2)

# Extract BUSCO
BUSCO_SUM=$(grep "C:" "${OUT}/qc/busco/short_summary.specific.bacteria_odb10.busco.txt" 2>/dev/null \
    | head -1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || echo "N/A")

# Extract Prokka
CDS=$(grep "^CDS" "${OUT}/prokka/MRSA.txt" | awk '{print $2}')
TRNA=$(grep "^tRNA" "${OUT}/prokka/MRSA.txt" | awk '{print $2}')
RRNA=$(grep "^rRNA" "${OUT}/prokka/MRSA.txt" | awk '{print $2}')

# IS elements
IS_COUNT=$(tail -n +2 "${OUT}/isescan/MRSA.fa.is.csv" 2>/dev/null | wc -l || echo "0")
IS_COUNT=$(echo $IS_COUNT | tr -d ' ')

# Integrons
INTEGRON_COUNT=$(grep -c "^" "${OUT}/integrons/Results_Integron_Finder_MRSA/MRSA.integrons" 2>/dev/null || echo "0")
INTEGRON_COUNT=$(echo $INTEGRON_COUNT | tr -d ' ')

# AMR counts
AMR_PROTEIN=$(tail -n +2 "${OUT}/amr/amrfinder_results.tsv" 2>/dev/null | wc -l || echo "0")
AMR_PROTEIN=$(echo $AMR_PROTEIN | tr -d ' ')
AMR_NUC=$(tail -n +2 "${OUT}/amr/abricate_card.tsv" 2>/dev/null | wc -l || echo "0")
AMR_NUC=$(echo $AMR_NUC | tr -d ' ')

# Virulence
VF_COUNT=$(tail -n +2 "${OUT}/virulence/abricate_vfdb.tsv" 2>/dev/null | wc -l || echo "0")
VF_COUNT=$(echo $VF_COUNT | tr -d ' ')

# Plasmids
PLASMID_COUNT=$(grep -c "plasmid" "${OUT}/plasmids/contig_report.txt" 2>/dev/null || echo "0")
PLASMID_COUNT=$(echo $PLASMID_COUNT | tr -d ' ')
REPLICON_COUNT=$(tail -n +2 "${OUT}/plasmids/plasmidfinder_results.tsv" 2>/dev/null | wc -l || echo "0")
REPLICON_COUNT=$(echo $REPLICON_COUNT | tr -d ' ')

# Write main CSV
cat > "${RES}/genomic_characterization.csv" << CSVEOF
metric,value
total_length,${TOTAL_LEN}
num_contigs,${NUM_CONTIGS}
n50,${N50}
gc_content,${GC}
largest_contig,${LARGEST}
completeness,${BUSCO_SUM}
cds_count,${CDS}
trna_count,${TRNA}
rrna_count,${RRNA}
is_elements_found,${IS_COUNT}
integrons_found,${INTEGRON_COUNT}
amr_genes_protein_method,${AMR_PROTEIN}
amr_genes_nucleotide_method,${AMR_NUC}
amr_genes_confirmed_both,${AMR_BOTH}
virulence_factors,${VF_COUNT}
plasmid_contigs,${PLASMID_COUNT}
replicon_types,${REPLICON_COUNT}
CSVEOF

# Write detailed AMR cross-validation CSV
cp "${OUT}/amr/cross_validated.csv" "${RES}/amr_cross_validated.csv"

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/genomic_characterization.csv"
echo ""
ls -lh "${RES}/"
