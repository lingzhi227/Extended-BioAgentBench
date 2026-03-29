#!/usr/bin/env bash
set -euo pipefail
#
# Circular RNA Discovery Pipeline
#
# DAG Structure (depth=8, convergence=3):
#
#   FASTQ_R1+R2 → trim(1) → STAR_index(2) → STAR_align(3)
#     ├→ CIRCexplorer2(4a) ─┐
#     └→ DCC(4b) ────────────┤ CONVERGE1 (consensus ≥2 tools)
#                             ↓
#                     consensus(5)
#       ├→ quantify(6a) ──────────┐
#       └→ extract_seqs(6b) ──────┤ CONVERGE2
#                                  ↓
#                          merge(7) → report(8)
#
# Inputs:  data/*_1.fastq.gz + *_2.fastq.gz, reference/chrI.fa, reference/chrI.gtf, reference/chrI.txt
# Outputs: results/report.csv

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
REF_DIR="${SCRIPT_DIR}/reference"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

GENOME="${REF_DIR}/chrI.fa"
GTF="${REF_DIR}/chrI.gtf"
GENEPRED="${REF_DIR}/chrI.txt"
STAR_IDX="${OUTPUT_DIR}/star_index"

mkdir -p "${OUTPUT_DIR}"/{trimmed,star_index,aligned,circexplorer2,dcc,consensus} "${RESULTS_DIR}"

echo "=== Circular RNA Discovery Pipeline ==="
echo "Threads: ${THREADS}"

# ============================================================
# Step 1: Build STAR index
# ============================================================
if [ ! -f "${STAR_IDX}/SA" ]; then
    echo "[Step 1] Building STAR genome index..."
    STAR --runMode genomeGenerate \
        --genomeDir "${STAR_IDX}" \
        --genomeFastaFiles "${GENOME}" \
        --sjdbGTFfile "${GTF}" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 11 \
        --runThreadN ${THREADS} 2>&1 | tail -5
else
    echo "[Step 1] STAR index exists, skipping."
fi

# ============================================================
# Step 2: Trim adapters for each sample
# ============================================================
SAMPLES=""
for r1 in "${DATA_DIR}"/*_1.fastq.gz; do
    SAMPLE=$(basename "$r1" _1.fastq.gz)
    SAMPLES="${SAMPLES} ${SAMPLE}"
    r2="${DATA_DIR}/${SAMPLE}_2.fastq.gz"
    if [ ! -f "${OUTPUT_DIR}/trimmed/${SAMPLE}_1_val_1.fq.gz" ]; then
        echo "[Step 2] Trimming: ${SAMPLE}..."
        trim_galore --paired --cores 4 \
            -o "${OUTPUT_DIR}/trimmed" \
            "$r1" "$r2" 2>&1 | tail -3
    fi
done

# ============================================================
# Step 3: STAR alignment with chimeric output
# ============================================================
for SAMPLE in ${SAMPLES}; do
    ALIGNED_DIR="${OUTPUT_DIR}/aligned/${SAMPLE}"
    if [ ! -f "${ALIGNED_DIR}/Chimeric.out.junction" ]; then
        echo "[Step 3] STAR alignment: ${SAMPLE}..."
        mkdir -p "${ALIGNED_DIR}"
        STAR --runThreadN ${THREADS} \
            --genomeDir "${STAR_IDX}" \
            --readFilesIn \
                "${OUTPUT_DIR}/trimmed/${SAMPLE}_1_val_1.fq.gz" \
                "${OUTPUT_DIR}/trimmed/${SAMPLE}_2_val_2.fq.gz" \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix "${ALIGNED_DIR}/" \
            --chimOutType Junctions SeparateSAMold \
            --chimJunctionOverhangMin 10 \
            --alignSJDBoverhangMin 10 \
            --chimSegmentMin 10 \
            --outSJtype Standard \
            --outSAMunmapped Within \
            --limitBAMsortRAM 2000000000 2>&1 | tail -5
        samtools index "${ALIGNED_DIR}/Aligned.sortedByCoord.out.bam"
    fi
done

# ============================================================
# Step 4a: CIRCexplorer2 detection
# ============================================================
for SAMPLE in ${SAMPLES}; do
    CE2_OUT="${OUTPUT_DIR}/circexplorer2/${SAMPLE}"
    if [ ! -f "${CE2_OUT}/circularRNA_known.txt" ]; then
        echo "[Step 4a] CIRCexplorer2: ${SAMPLE}..."
        mkdir -p "${CE2_OUT}"
        # Parse STAR chimeric junctions
        CIRCexplorer2 parse -t STAR \
            "${OUTPUT_DIR}/aligned/${SAMPLE}/Chimeric.out.junction" \
            -b "${CE2_OUT}/back_spliced_junction.bed" 2>&1 | tail -3
        # Annotate with gene model
        CIRCexplorer2 annotate \
            -r "${GENEPRED}" \
            -g "${GENOME}" \
            -b "${CE2_OUT}/back_spliced_junction.bed" \
            -o "${CE2_OUT}/circularRNA_known.txt" 2>&1 | tail -3
    fi
done

# ============================================================
# Step 4b: DCC detection
# ============================================================
DCC_OUT="${OUTPUT_DIR}/dcc"
if [ ! -f "${DCC_OUT}/CircRNACount" ]; then
    echo "[Step 4b] DCC detection..."
    mkdir -p "${DCC_OUT}"
    # Prepare file lists for DCC
    MATE1_LIST="${DCC_OUT}/samplesheet_mate1.txt"
    MATE2_LIST="${DCC_OUT}/samplesheet_mate2.txt"
    CHIMERIC_LIST="${DCC_OUT}/samplesheet_chimeric.txt"
    > "${MATE1_LIST}"
    > "${MATE2_LIST}"
    > "${CHIMERIC_LIST}"
    for SAMPLE in ${SAMPLES}; do
        echo "${OUTPUT_DIR}/aligned/${SAMPLE}/Chimeric.out.junction" >> "${CHIMERIC_LIST}"
    done

    cd "${DCC_OUT}"
    DCC "${CHIMERIC_LIST}" \
        -an "${GTF}" \
        -F -M -Nr 1 1 \
        -fg -G -A "${GENOME}" \
        -D -T ${THREADS} 2>&1 | tail -10 || true
    cd "${SCRIPT_DIR}"
else
    echo "[Step 4b] DCC already done, skipping."
fi

# ============================================================
# Step 5: Consensus (CONVERGE1 — ≥2 tools must agree)
# ============================================================
echo "[Step 5] Building consensus..."
python3 "${SCRIPT_DIR}/scripts/consensus.py" \
    "${OUTPUT_DIR}/circexplorer2" \
    "${OUTPUT_DIR}/dcc" \
    "${OUTPUT_DIR}/consensus"

# ============================================================
# Step 6a+6b: Quantify + extract sequences (CONVERGE2)
# ============================================================
echo "[Step 6] Quantifying and extracting sequences..."
python3 "${SCRIPT_DIR}/scripts/quantify_and_report.py" \
    "${OUTPUT_DIR}/consensus/consensus_circrnas.bed" \
    "${OUTPUT_DIR}/circexplorer2" \
    "${GENOME}" \
    "${RESULTS_DIR}"

# ============================================================
# Validation
# ============================================================
echo ""
echo "=== Validation ==="
if [ -f "${RESULTS_DIR}/report.csv" ]; then
    N_METRICS=$(tail -n +2 "${RESULTS_DIR}/report.csv" | wc -l)
    N_EMPTY=$(grep -cE ",,|,$" "${RESULTS_DIR}/report.csv" || true)
    echo "Report: ${N_METRICS} metrics, ${N_EMPTY} empty values"
    echo ""
    cat "${RESULTS_DIR}/report.csv"
else
    echo "ERROR: Report not generated!"
    exit 1
fi

echo ""
echo "=== Pipeline complete ==="
