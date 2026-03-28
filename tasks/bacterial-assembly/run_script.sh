#!/bin/bash
set -e

# =============================================================================
# Bacterial Genome Assembly + Multi-Tool Annotation (Diamond DAG)
#
# DAG structure:
#   FASTQ → fastp → shovill (assemble)
#                        ↓
#           ┌────────┬───┼────┬─────────┐
#         Prokka  ABRicate  mlst  BUSCO  QUAST
#         (genes) (AMR)   (type) (QC)   (stats)
#           └────────┴───┼────┴─────────┘
#                        ↓
#                  Merge → CSV report
#
# Data: MRSA (S. aureus) Illumina paired-end (Hikichi et al. 2019)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "=================================================================="
    echo "Started at: $(date)"
}

mkdir -p "${OUT_DIR}/trimmed" "${OUT_DIR}/assembly" "${OUT_DIR}/prokka" \
         "${OUT_DIR}/abricate" "${OUT_DIR}/mlst" "${OUT_DIR}/busco" "${OUT_DIR}/quast"
mkdir -p "${RESULTS_DIR}"

# ==========================================================================
# STEP 1: Quality trimming with fastp
# ==========================================================================
log_step "Quality trimming with fastp"

if [ ! -f "${OUT_DIR}/trimmed/R1.fastq" ]; then
    fastp \
        --in1 "${DATA_DIR}/reads_R1.fastq" \
        --in2 "${DATA_DIR}/reads_R2.fastq" \
        --out1 "${OUT_DIR}/trimmed/R1.fastq" \
        --out2 "${OUT_DIR}/trimmed/R2.fastq" \
        --detect_adapter_for_pe \
        --cut_front --cut_tail \
        --cut_mean_quality 20 \
        --length_required 30 \
        --thread ${THREADS} \
        --json "${OUT_DIR}/trimmed/fastp.json" \
        --html "${OUT_DIR}/trimmed/fastp.html"
    echo "Trimming complete."
else
    echo "Trimmed reads exist, skipping."
fi

# ==========================================================================
# STEP 2: Genome assembly with shovill (SPAdes wrapper)
# ==========================================================================
log_step "Genome assembly with shovill"

if [ ! -f "${OUT_DIR}/assembly/contigs.fa" ]; then
    shovill \
        --R1 "${OUT_DIR}/trimmed/R1.fastq" \
        --R2 "${OUT_DIR}/trimmed/R2.fastq" \
        --outdir "${OUT_DIR}/assembly" \
        --gsize 2914567 \
        --cpus ${THREADS} \
        --ram 8 \
        --force
    echo "Assembly complete."
else
    echo "Assembly exists, skipping."
fi

CONTIGS="${OUT_DIR}/assembly/contigs.fa"

# === DIAMOND BRANCHES START HERE ===
# The following 5 tools run independently on the same assembly

# ==========================================================================
# BRANCH 1: Gene annotation with Prokka
# ==========================================================================
log_step "BRANCH 1: Gene annotation with Prokka"

if [ ! -f "${OUT_DIR}/prokka/MRSA.gff" ]; then
    prokka \
        "${CONTIGS}" \
        --outdir "${OUT_DIR}/prokka" \
        --prefix MRSA \
        --cpus ${THREADS} \
        --kingdom Bacteria \
        --genus Staphylococcus \
        --species aureus \
        --force
    echo "Prokka complete."
else
    echo "Prokka output exists, skipping."
fi

# ==========================================================================
# BRANCH 2: AMR gene detection with ABRicate
# ==========================================================================
log_step "BRANCH 2: AMR gene detection with ABRicate"

if [ ! -f "${OUT_DIR}/abricate/card_results.tsv" ]; then
    # CARD database for AMR
    abricate "${CONTIGS}" --db card --minid 80 --mincov 60 \
        > "${OUT_DIR}/abricate/card_results.tsv"
    # VFDB for virulence factors
    abricate "${CONTIGS}" --db vfdb --minid 80 --mincov 60 \
        > "${OUT_DIR}/abricate/vfdb_results.tsv"
    echo "ABRicate complete."
else
    echo "ABRicate output exists, skipping."
fi

# ==========================================================================
# BRANCH 3: MLST typing
# ==========================================================================
log_step "BRANCH 3: MLST typing"

if [ ! -f "${OUT_DIR}/mlst/mlst_results.tsv" ]; then
    mlst "${CONTIGS}" > "${OUT_DIR}/mlst/mlst_results.tsv"
    echo "MLST complete."
    cat "${OUT_DIR}/mlst/mlst_results.tsv"
else
    echo "MLST output exists, skipping."
fi

# ==========================================================================
# BRANCH 4: Assembly completeness with BUSCO
# ==========================================================================
log_step "BRANCH 4: Assembly completeness with BUSCO"

if [ ! -f "${OUT_DIR}/busco/short_summary.specific.bacteria_odb10.busco.txt" ]; then
    busco \
        -i "${CONTIGS}" \
        -o busco \
        --out_path "${OUT_DIR}" \
        -l bacteria_odb10 \
        -m genome \
        -c ${THREADS} \
        --force
    echo "BUSCO complete."
else
    echo "BUSCO output exists, skipping."
fi

# ==========================================================================
# BRANCH 5: Assembly statistics with QUAST
# ==========================================================================
log_step "BRANCH 5: Assembly statistics with QUAST"

if [ ! -f "${OUT_DIR}/quast/report.tsv" ]; then
    quast "${CONTIGS}" -o "${OUT_DIR}/quast" -t ${THREADS}
    echo "QUAST complete."
else
    echo "QUAST output exists, skipping."
fi

# === DIAMOND MERGE ===

# ==========================================================================
# STEP 8: Merge all branch results into comprehensive CSV
# ==========================================================================
log_step "MERGE: Generating comprehensive results CSV"

# --- Extract QUAST metrics ---
TOTAL_LENGTH=$(grep "^Total length" "${OUT_DIR}/quast/report.tsv" | head -1 | cut -f2)
NUM_CONTIGS=$(grep "^# contigs " "${OUT_DIR}/quast/report.tsv" | head -1 | cut -f2)
N50=$(grep "^N50" "${OUT_DIR}/quast/report.tsv" | cut -f2)
GC_CONTENT=$(grep "^GC" "${OUT_DIR}/quast/report.tsv" | cut -f2)
LARGEST_CONTIG=$(grep "^Largest contig" "${OUT_DIR}/quast/report.tsv" | cut -f2)

# --- Extract BUSCO summary ---
BUSCO_SUMMARY=$(grep "C:" "${OUT_DIR}/busco/short_summary.specific.bacteria_odb10.busco.txt" 2>/dev/null \
    | head -1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' \
    || echo "N/A")

# --- Extract MLST ---
MLST_SCHEME=$(cut -f2 "${OUT_DIR}/mlst/mlst_results.tsv")
MLST_ST=$(cut -f3 "${OUT_DIR}/mlst/mlst_results.tsv")

# --- Extract Prokka gene counts ---
PROKKA_CDS=$(grep "^CDS" "${OUT_DIR}/prokka/MRSA.txt" 2>/dev/null | awk '{print $2}' || echo "0")
PROKKA_TRNA=$(grep "^tRNA" "${OUT_DIR}/prokka/MRSA.txt" 2>/dev/null | awk '{print $2}' || echo "0")
PROKKA_RRNA=$(grep "^rRNA" "${OUT_DIR}/prokka/MRSA.txt" 2>/dev/null | awk '{print $2}' || echo "0")

# --- Write assembly_report.csv (main output) ---
echo "metric,value" > "${RESULTS_DIR}/assembly_report.csv"
echo "total_length,${TOTAL_LENGTH}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "num_contigs,${NUM_CONTIGS}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "n50,${N50}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "gc_content,${GC_CONTENT}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "largest_contig,${LARGEST_CONTIG}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "busco_summary,${BUSCO_SUMMARY}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "mlst_scheme,${MLST_SCHEME}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "mlst_st,${MLST_ST}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "prokka_cds,${PROKKA_CDS}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "prokka_trna,${PROKKA_TRNA}" >> "${RESULTS_DIR}/assembly_report.csv"
echo "prokka_rrna,${PROKKA_RRNA}" >> "${RESULTS_DIR}/assembly_report.csv"

# --- Write AMR results CSV ---
echo "gene,accession,database,identity,coverage,contig" > "${RESULTS_DIR}/amr_genes.csv"
if [ -s "${OUT_DIR}/abricate/card_results.tsv" ]; then
    tail -n +2 "${OUT_DIR}/abricate/card_results.tsv" | \
        awk -F'\t' 'BEGIN{OFS=","} {print $6,$12,"CARD",$10,$11,$2}' \
        >> "${RESULTS_DIR}/amr_genes.csv"
fi

# --- Write virulence results CSV ---
echo "gene,accession,database,identity,coverage,contig" > "${RESULTS_DIR}/virulence_genes.csv"
if [ -s "${OUT_DIR}/abricate/vfdb_results.tsv" ]; then
    tail -n +2 "${OUT_DIR}/abricate/vfdb_results.tsv" | \
        awk -F'\t' 'BEGIN{OFS=","} {print $6,$12,"VFDB",$10,$11,$2}' \
        >> "${RESULTS_DIR}/virulence_genes.csv"
fi

echo ""
echo "=== Pipeline complete ==="
echo "=== assembly_report.csv ==="
cat "${RESULTS_DIR}/assembly_report.csv"
echo ""
echo "=== AMR genes found ==="
wc -l "${RESULTS_DIR}/amr_genes.csv"
head -5 "${RESULTS_DIR}/amr_genes.csv"
echo ""
echo "=== Virulence genes found ==="
wc -l "${RESULTS_DIR}/virulence_genes.csv"
head -5 "${RESULTS_DIR}/virulence_genes.csv"
echo ""
ls -lh "${RESULTS_DIR}/"
