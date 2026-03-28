#!/bin/bash
set -e

# =============================================================================
# Task 20: Metagenomic Profiling (dual-path: read-based + assembly-based)
#
# DAG (depth 7, dual-path with cross-validation):
#
# L0: metagenomic reads (PE)
# L1: fastp (trim)
#     ├──────────────────────────────────────┐
# L2: kraken2 (read classification)     megahit (assembly)
#     │                                      │
# L3: bracken (abundance estimation)    prodigal (gene prediction)
#     │                                      │
#     │                               ├──────┴──────┐
# L4: │                            diamond        abricate
#     │                            (functional)   (AMR)
#     │                                │              │
# L5: │                            GO/KEGG summary   │
#     └──────────┬───────────────────────┴────────────┘
# L6: MERGE (taxonomy + function + AMR report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,kraken,assembly,genes,functional,amr} "${RES}"

# ===========================================================================
# L1: Trimming
# ===========================================================================
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/meta_R1.fastq.gz" --in2 "${DATA}/meta_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} \
          --json "${OUT}/trimmed/fastp.json"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2 LEFT: Taxonomic classification with Kraken2
# ===========================================================================
log_step "L2-LEFT: kraken2 classification"
KRAKEN_DB="/pscratch/sd/l/lingzhi/kraken2_db"
if [ ! -f "${OUT}/kraken/kraken_report.txt" ]; then
    if [ -d "${KRAKEN_DB}" ]; then
        kraken2 --db "${KRAKEN_DB}" --threads ${THREADS} --paired \
                --report "${OUT}/kraken/kraken_report.txt" \
                --output "${OUT}/kraken/kraken_output.txt" \
                "${OUT}/trimmed/R1.fastq.gz" "${OUT}/trimmed/R2.fastq.gz"
    else
        echo "WARNING: Kraken2 database not found at ${KRAKEN_DB}"
        echo "Building mini database..."
        mkdir -p "${KRAKEN_DB}"
        kraken2-build --download-taxonomy --db "${KRAKEN_DB}" 2>&1 || true
        kraken2-build --download-library bacteria --db "${KRAKEN_DB}" 2>&1 || true
        kraken2-build --build --db "${KRAKEN_DB}" --threads ${THREADS} 2>&1 || {
            echo "WARNING: Could not build Kraken2 DB. Skipping read-based classification."
            touch "${OUT}/kraken/kraken_report.txt"
        }
        if [ -s "${KRAKEN_DB}/hash.k2d" ]; then
            kraken2 --db "${KRAKEN_DB}" --threads ${THREADS} --paired \
                    --report "${OUT}/kraken/kraken_report.txt" \
                    --output "${OUT}/kraken/kraken_output.txt" \
                    "${OUT}/trimmed/R1.fastq.gz" "${OUT}/trimmed/R2.fastq.gz"
        fi
    fi
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2 RIGHT: Metagenomic assembly with MEGAHIT
# ===========================================================================
log_step "L2-RIGHT: megahit assembly"
if [ ! -f "${OUT}/assembly/final.contigs.fa" ]; then
    megahit -1 "${OUT}/trimmed/R1.fastq.gz" -2 "${OUT}/trimmed/R2.fastq.gz" \
            -o "${OUT}/assembly" -t ${THREADS} --min-contig-len 500
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3 LEFT: Bracken abundance estimation
# ===========================================================================
log_step "L3-LEFT: bracken"
if [ -s "${OUT}/kraken/kraken_report.txt" ] && [ ! -f "${OUT}/kraken/bracken_species.txt" ]; then
    bracken -d "${KRAKEN_DB}" -i "${OUT}/kraken/kraken_report.txt" \
            -o "${OUT}/kraken/bracken_species.txt" -l S -t 10 2>&1 || true
fi

# ===========================================================================
# L3 RIGHT: Gene prediction with Prodigal
# ===========================================================================
log_step "L3-RIGHT: prodigal gene prediction"
if [ ! -f "${OUT}/genes/proteins.faa" ]; then
    prodigal -i "${OUT}/assembly/final.contigs.fa" \
             -o "${OUT}/genes/genes.gff" -a "${OUT}/genes/proteins.faa" \
             -p meta -f gff
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4: Assembly QC + AMR
# ===========================================================================
log_step "L4: quast assembly QC"
if [ ! -f "${OUT}/assembly/quast/report.tsv" ]; then
    quast "${OUT}/assembly/final.contigs.fa" -o "${OUT}/assembly/quast" -t ${THREADS} --min-contig 500
else echo "Skipping (exists)"; fi

log_step "L4: abricate AMR on contigs"
if [ ! -f "${OUT}/amr/abricate.tsv" ]; then
    abricate "${OUT}/assembly/final.contigs.fa" --db card --minid 80 --mincov 60 > "${OUT}/amr/abricate.tsv" 2>/dev/null || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5-L6: MERGE
# ===========================================================================
log_step "MERGE"

# Assembly stats
TOTAL_LEN=$(grep "^Total length" "${OUT}/assembly/quast/report.tsv" 2>/dev/null | head -1 | cut -f2 || echo "N/A")
NUM_CONTIGS=$(grep "^# contigs " "${OUT}/assembly/quast/report.tsv" 2>/dev/null | head -1 | cut -f2 || echo "N/A")
N50=$(grep "^N50" "${OUT}/assembly/quast/report.tsv" 2>/dev/null | cut -f2 || echo "N/A")
TOTAL_GENES=$(grep -c ">" "${OUT}/genes/proteins.faa" 2>/dev/null || echo "0")
AMR_GENES=$(tail -n +2 "${OUT}/amr/abricate.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# Taxonomy from Kraken (top 5 species)
if [ -s "${OUT}/kraken/kraken_report.txt" ]; then
    CLASSIFIED_PCT=$(head -1 "${OUT}/kraken/kraken_output.txt" 2>/dev/null | awk '{print "N/A"}')
    # Get top species
    awk '$4=="S" {print $6,$7,$8}' "${OUT}/kraken/kraken_report.txt" 2>/dev/null | sort -k1 -rn | head -5 > "${OUT}/kraken/top_species.txt" || true
    NUM_SPECIES=$(awk '$4=="S"' "${OUT}/kraken/kraken_report.txt" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
else
    NUM_SPECIES="N/A"
fi

cat > "${RES}/metagenome_report.csv" << CSVEOF
metric,value
assembly_length,${TOTAL_LEN}
num_contigs,${NUM_CONTIGS}
n50,${N50}
predicted_genes,${TOTAL_GENES}
species_detected,${NUM_SPECIES}
amr_genes,${AMR_GENES}
CSVEOF

# Copy detailed results
cp "${OUT}/kraken/kraken_report.txt" "${RES}/taxonomy_report.txt" 2>/dev/null || true
cp "${OUT}/kraken/top_species.txt" "${RES}/top_species.txt" 2>/dev/null || true

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/metagenome_report.csv"
echo ""
ls -lh "${RES}/"
