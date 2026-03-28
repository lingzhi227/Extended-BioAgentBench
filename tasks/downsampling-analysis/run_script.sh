#!/bin/bash
set -e

# =============================================================================
# Task 27: Read Downsampling and Assembly Quality Titration
#
# DAG (depth 5, replicated fan-out):
#
# L0: PE reads + reference
# L1: fastp (trim full dataset)
# L2: ├── seqkit sample 25% → megahit → quast
#     ├── seqkit sample 50% → megahit → quast
#     └── seqkit sample 100% → megahit → quast
# L3: Parse all three QUAST reports
# L4: MERGE (titration report: how coverage affects assembly)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

REFERENCE="${REF}/reference.fna"
FRACTIONS=("0.25" "0.50" "1.00")

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}/trimmed" "${RES}"

# L1: Trim
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} --json "${OUT}/trimmed/fastp.json"
else echo "Skipping (exists)"; fi

# L2: For each fraction, subsample + assemble + QC
for FRAC in "${FRACTIONS[@]}"; do
    FRAC_DIR="${OUT}/frac_${FRAC}"
    mkdir -p "${FRAC_DIR}"

    log_step "L2: Subsample ${FRAC} + assemble"
    if [ ! -f "${FRAC_DIR}/assembly/final.contigs.fa" ]; then
        if [ "${FRAC}" = "1.00" ]; then
            cp "${OUT}/trimmed/R1.fastq.gz" "${FRAC_DIR}/R1.fastq.gz"
            cp "${OUT}/trimmed/R2.fastq.gz" "${FRAC_DIR}/R2.fastq.gz"
        else
            seqkit sample -p "${FRAC}" -s 42 "${OUT}/trimmed/R1.fastq.gz" -o "${FRAC_DIR}/R1.fastq.gz"
            seqkit sample -p "${FRAC}" -s 42 "${OUT}/trimmed/R2.fastq.gz" -o "${FRAC_DIR}/R2.fastq.gz"
        fi
        megahit -1 "${FRAC_DIR}/R1.fastq.gz" -2 "${FRAC_DIR}/R2.fastq.gz" \
                -o "${FRAC_DIR}/assembly" -t ${THREADS} --min-contig-len 500
    else echo "Skipping (exists)"; fi

    log_step "L2: QUAST ${FRAC}"
    if [ ! -f "${FRAC_DIR}/quast/report.tsv" ]; then
        quast "${FRAC_DIR}/assembly/final.contigs.fa" -r "${REFERENCE}" \
              -o "${FRAC_DIR}/quast" -t ${THREADS}
    else echo "Skipping (exists)"; fi
done

# L3-L4: Parse + MERGE
log_step "MERGE"

echo "fraction,total_length,num_contigs,n50,largest_contig,genome_fraction,num_reads" > "${RES}/downsampling_report.csv"

for FRAC in "${FRACTIONS[@]}"; do
    FRAC_DIR="${OUT}/frac_${FRAC}"
    TOTAL=$(grep "^Total length" "${FRAC_DIR}/quast/report.tsv" | head -1 | cut -f2)
    NCTG=$(grep "^# contigs " "${FRAC_DIR}/quast/report.tsv" | head -1 | cut -f2)
    N50=$(grep "^N50" "${FRAC_DIR}/quast/report.tsv" | cut -f2)
    LARGEST=$(grep "^Largest contig" "${FRAC_DIR}/quast/report.tsv" | cut -f2)
    GF=$(grep "^Genome fraction" "${FRAC_DIR}/quast/report.tsv" | cut -f2)
    NREADS=$(zcat "${FRAC_DIR}/R1.fastq.gz" | wc -l | awk '{print $1/4}')
    echo "${FRAC},${TOTAL},${NCTG},${N50},${LARGEST},${GF},${NREADS}" >> "${RES}/downsampling_report.csv"
done

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/downsampling_report.csv"
