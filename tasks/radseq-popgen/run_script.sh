#!/usr/bin/env bash
set -euo pipefail
#
# RADseq Population Genetics Pipeline (de novo)
#
# DAG Structure (depth=9, convergence=2):
#
#   demultiplex(1)
#     → per_sample_stacks(2) [parallel per sample]
#       → build_catalog(3)  ← CONVERGE1 (all sample stacks → catalog)
#         → match_catalog(4)
#           → transpose_data(5)
#             → genotype_calling(6)
#               → population_stats(7) ← CONVERGE2 (genotypes + pop map)
#                 ├→ export_vcf(8a)
#                 └→ export_structure(8b)
#                   → report(9)
#
# Inputs:  data/SRR034310.fastq.gz, reference/barcodes.txt, reference/popmap.txt
# Outputs: results/report.csv

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
REF_DIR="${SCRIPT_DIR}/reference"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

DEMUX_DIR="${OUTPUT_DIR}/demultiplexed"
STACKS_DIR="${OUTPUT_DIR}/stacks"

mkdir -p "${DEMUX_DIR}" "${STACKS_DIR}" "${RESULTS_DIR}"

echo "=== RADseq Population Genetics Pipeline ==="
echo "Threads: ${THREADS}"

# ============================================================
# Step 1: Demultiplex raw reads using barcodes
# ============================================================
if [ ! -f "${DEMUX_DIR}/.done" ]; then
    echo "[Step 1] Demultiplexing reads..."
    process_radtags \
        -f "${DATA_DIR}/SRR034310.fastq.gz" \
        -o "${DEMUX_DIR}" \
        -b "${REF_DIR}/barcodes.txt" \
        --inline_null \
        -e sbfI \
        -r -c -q \
        --renz_1 sbfI \
        -i gzfastq 2>&1 | tee "${OUTPUT_DIR}/demux.log"
    touch "${DEMUX_DIR}/.done"
else
    echo "[Step 1] Demultiplexing already done, skipping."
fi

# Count demultiplexed samples
N_SAMPLES=$(ls "${DEMUX_DIR}"/*.fq.gz 2>/dev/null | grep -v remainder | wc -l)
echo "Demultiplexed samples: ${N_SAMPLES}"

if [ "${N_SAMPLES}" -lt 2 ]; then
    echo "ERROR: Too few samples demultiplexed (${N_SAMPLES}). Check barcodes."
    exit 1
fi

# ============================================================
# Step 2: Build per-sample loci (ustacks)
# ============================================================
if [ ! -f "${STACKS_DIR}/.ustacks_done" ]; then
    echo "[Step 2] Building per-sample loci..."
    SAMPLE_ID=1
    for fq in "${DEMUX_DIR}"/*.fq.gz; do
        BASENAME=$(basename "$fq" .fq.gz)
        # Skip remainder file
        if [[ "$BASENAME" == *"remainder"* ]]; then continue; fi
        echo "  ustacks: ${BASENAME} (id=${SAMPLE_ID})"
        ustacks -f "$fq" -o "${STACKS_DIR}" -i ${SAMPLE_ID} \
            -M 3 -m 3 -p ${THREADS} 2>&1 | tail -2
        SAMPLE_ID=$((SAMPLE_ID + 1))
    done
    touch "${STACKS_DIR}/.ustacks_done"
else
    echo "[Step 2] ustacks already done, skipping."
fi

# ============================================================
# Step 3: Build catalog (CONVERGE1 — all sample stacks merge)
# ============================================================
if [ ! -f "${STACKS_DIR}/.cstacks_done" ]; then
    echo "[Step 3] Building catalog..."
    # Build -s flags for all samples
    SAMPLE_FLAGS=""
    for fq in "${DEMUX_DIR}"/*.fq.gz; do
        BASENAME=$(basename "$fq" .fq.gz)
        if [[ "$BASENAME" == *"remainder"* ]]; then continue; fi
        SAMPLE_FLAGS="${SAMPLE_FLAGS} -s ${STACKS_DIR}/${BASENAME}"
    done
    cstacks ${SAMPLE_FLAGS} -o "${STACKS_DIR}" -n 3 -p ${THREADS} 2>&1 | tail -5
    touch "${STACKS_DIR}/.cstacks_done"
else
    echo "[Step 3] Catalog already built, skipping."
fi

# ============================================================
# Step 4: Match samples against catalog (sstacks)
# ============================================================
if [ ! -f "${STACKS_DIR}/.sstacks_done" ]; then
    echo "[Step 4] Matching samples against catalog..."
    for fq in "${DEMUX_DIR}"/*.fq.gz; do
        BASENAME=$(basename "$fq" .fq.gz)
        if [[ "$BASENAME" == *"remainder"* ]]; then continue; fi
        echo "  sstacks: ${BASENAME}"
        sstacks -c "${STACKS_DIR}" -s "${STACKS_DIR}/${BASENAME}" \
            -o "${STACKS_DIR}" -p ${THREADS} 2>&1 | tail -2
    done
    touch "${STACKS_DIR}/.sstacks_done"
else
    echo "[Step 4] sstacks already done, skipping."
fi

# ============================================================
# Step 5: Transpose data (tsv2bam)
# ============================================================
if [ ! -f "${STACKS_DIR}/.tsv2bam_done" ]; then
    echo "[Step 5] Transposing data..."
    tsv2bam -P "${STACKS_DIR}" -M "${REF_DIR}/popmap.txt" -t ${THREADS} 2>&1 | tail -5
    touch "${STACKS_DIR}/.tsv2bam_done"
else
    echo "[Step 5] tsv2bam already done, skipping."
fi

# ============================================================
# Step 6: Genotype calling (gstacks)
# ============================================================
if [ ! -f "${STACKS_DIR}/.gstacks_done" ]; then
    echo "[Step 6] Calling genotypes..."
    gstacks -P "${STACKS_DIR}" -M "${REF_DIR}/popmap.txt" -t ${THREADS} 2>&1 | tail -5
    touch "${STACKS_DIR}/.gstacks_done"
else
    echo "[Step 6] gstacks already done, skipping."
fi

# ============================================================
# Step 7: Population statistics (CONVERGE2 — genotypes + pop map)
# Step 8a: Export VCF
# Step 8b: Export STRUCTURE
# ============================================================
if [ ! -f "${OUTPUT_DIR}/populations/.done" ]; then
    echo "[Step 7-8] Computing population statistics..."
    mkdir -p "${OUTPUT_DIR}/populations"
    populations -P "${STACKS_DIR}" -M "${REF_DIR}/popmap.txt" \
        -O "${OUTPUT_DIR}/populations" \
        -p 2 -r 0.8 \
        --vcf --structure --fstats \
        -t ${THREADS} 2>&1 | tail -10
    touch "${OUTPUT_DIR}/populations/.done"
else
    echo "[Step 7-8] Populations already done, skipping."
fi

# ============================================================
# Step 9: Generate report
# ============================================================
echo "[Step 9] Generating report..."

# Parse demux log
TOTAL_READS=$(grep -oP 'total sequences\s+\K\d+' "${OUTPUT_DIR}/demux.log" 2>/dev/null || echo "0")
if [ "$TOTAL_READS" = "0" ]; then
    TOTAL_READS=$(grep -m1 "total sequences" "${OUTPUT_DIR}/demux.log" | grep -oP '\d+' | head -1 || echo "0")
fi
RETAINED_READS=$(grep -oP 'retained reads\s+\K\d+' "${OUTPUT_DIR}/demux.log" 2>/dev/null || echo "0")
if [ "$RETAINED_READS" = "0" ]; then
    RETAINED_READS=$(grep -m1 "retained" "${OUTPUT_DIR}/demux.log" | grep -oP '\d[\d,]+' | head -1 | tr -d ',' || echo "0")
fi

# Count catalog loci
CATALOG_LOCI=$(zcat "${STACKS_DIR}/catalog.tags.tsv.gz" 2>/dev/null | tail -n +2 | wc -l || echo "0")
if [ "$CATALOG_LOCI" = "0" ]; then
    CATALOG_LOCI=$(ls "${STACKS_DIR}"/catalog.* 2>/dev/null | head -1 | xargs -I{} sh -c 'zcat {} 2>/dev/null | wc -l' || echo "0")
fi

# Parse populations output
POP_LOG="${OUTPUT_DIR}/populations/populations.log"
VARIANT_SITES=0
LOCI_KEPT=0
if [ -f "$POP_LOG" ]; then
    VARIANT_SITES=$(grep -oP 'Kept \K\d+' "$POP_LOG" 2>/dev/null | tail -1 || echo "0")
    LOCI_KEPT=$(grep -oP 'Kept \K\d+' "$POP_LOG" 2>/dev/null | head -1 || echo "0")
fi

# Parse VCF for stats
VCF_FILE=$(ls "${OUTPUT_DIR}"/populations/*.vcf 2>/dev/null | head -1)
VCF_SNPS=0
VCF_SAMPLES=0
if [ -n "$VCF_FILE" ] && [ -f "$VCF_FILE" ]; then
    VCF_SNPS=$(grep -c -v "^#" "$VCF_FILE" || true)
    VCF_SAMPLES=$(grep "^#CHROM" "$VCF_FILE" | tr '\t' '\n' | tail -n +10 | wc -l || echo "0")
fi

# Parse STRUCTURE file
STRUCT_FILE=$(ls "${OUTPUT_DIR}"/populations/*.structure 2>/dev/null | head -1)
STRUCT_LOCI=0
if [ -n "$STRUCT_FILE" ] && [ -f "$STRUCT_FILE" ]; then
    STRUCT_LOCI=$(head -1 "$STRUCT_FILE" | tr '\t' '\n' | wc -l || echo "0")
fi

# Parse Fst
FST_FILE=$(ls "${OUTPUT_DIR}"/populations/*.fst_summary.tsv 2>/dev/null | head -1)
MEAN_FST="NA"
if [ -n "$FST_FILE" ] && [ -f "$FST_FILE" ]; then
    MEAN_FST=$(awk -F'\t' 'NR>1 {sum+=$NF; n++} END {if(n>0) printf "%.4f", sum/n; else print "NA"}' "$FST_FILE" 2>/dev/null || echo "NA")
fi

# Parse per-population nucleotide diversity
HAPSTATS_FILE=$(ls "${OUTPUT_DIR}"/populations/*.hapstats.tsv 2>/dev/null | head -1)
MEAN_PI="NA"
if [ -n "$HAPSTATS_FILE" ] && [ -f "$HAPSTATS_FILE" ]; then
    MEAN_PI=$(awk -F'\t' 'NR>1 && $NF != "-" {sum+=$NF; n++} END {if(n>0) printf "%.6f", sum/n; else print "NA"}' "$HAPSTATS_FILE" 2>/dev/null || echo "NA")
fi

# Write report
cat > "${RESULTS_DIR}/report.csv" << CSVEOF
metric,value
total_reads,${TOTAL_READS}
samples_demultiplexed,${N_SAMPLES}
catalog_loci,${CATALOG_LOCI}
variant_sites_vcf,${VCF_SNPS}
samples_in_vcf,${VCF_SAMPLES}
structure_loci,${STRUCT_LOCI}
loci_passing_filters,${LOCI_KEPT}
mean_fst,${MEAN_FST}
mean_nucleotide_diversity,${MEAN_PI}
populations,2
CSVEOF

echo ""
echo "=== Final Report ==="
cat "${RESULTS_DIR}/report.csv"

# Validate
N_EMPTY=$(grep -cE ",,|,$|,NA$|,0$" "${RESULTS_DIR}/report.csv" || true)
echo ""
echo "Empty/zero/NA values: ${N_EMPTY}"

echo ""
echo "=== Pipeline complete ==="
