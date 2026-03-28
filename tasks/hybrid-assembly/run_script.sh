#!/bin/bash
set -e

# =============================================================================
# Task 17: Hybrid Genome Assembly (Illumina + Nanopore)
#
# DAG (depth 6, V-merge — two input types converge):
#
# L0: illumina reads        nanopore reads
# L1: fastp (trim PE)       NanoPlot (QC) + Filtlong (filter)
#     └─────────┬───────────┘
# L2:      Unicycler (hybrid assembly using both)
#          ├──────────────────────────────┐
# L3:      │                      minimap2 (map nanopore → assembly)
#          │                         │
# L4:   QUAST (QC vs ref)      Polypolish (short-read polishing)
#        + BUSCO                     │
#          │                      Prokka (annotation)
#          │                         │
# L5:      │                      abricate + mlst
#          │                         │
# L6:   MERGE ───────────────────────┘
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

ILLUMINA_R1="${DATA}/illumina_R1.fastq"
ILLUMINA_R2="${DATA}/illumina_R2.fastq"
NANOPORE="${DATA}/nanopore.fastq"
REFERENCE="${REF}/reference.fasta"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,nanoqc,filtered,assembly,polished,qc,prokka,amr,mlst_out} "${RES}"

# ===========================================================================
# L1 LEFT: Illumina trimming
# ===========================================================================
log_step "L1-LEFT: fastp Illumina trimming"
if [ ! -f "${OUT}/trimmed/R1.fastq" ]; then
    fastp --in1 "${ILLUMINA_R1}" --in2 "${ILLUMINA_R2}" \
          --out1 "${OUT}/trimmed/R1.fastq" --out2 "${OUT}/trimmed/R2.fastq" \
          --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
          --length_required 50 --thread ${THREADS} \
          --json "${OUT}/trimmed/fastp.json" --html "${OUT}/trimmed/fastp.html"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L1 RIGHT: Nanopore QC + filtering
# ===========================================================================
log_step "L1-RIGHT: Nanopore QC + filtering"
if [ ! -f "${OUT}/filtered/nanopore_filtered.fastq" ]; then
    NanoPlot --fastq "${NANOPORE}" -o "${OUT}/nanoqc" -t ${THREADS} --plots dot 2>&1 || true
    filtlong --min_length 500 "${NANOPORE}" > "${OUT}/filtered/nanopore_filtered.fastq"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2: Hybrid assembly with Unicycler (V-merge: both data types converge)
# ===========================================================================
log_step "L2: Unicycler hybrid assembly"
if [ ! -f "${OUT}/assembly/assembly.fasta" ]; then
    unicycler -1 "${OUT}/trimmed/R1.fastq" -2 "${OUT}/trimmed/R2.fastq" \
              -l "${OUT}/filtered/nanopore_filtered.fastq" \
              -o "${OUT}/assembly" -t ${THREADS} --mode normal
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3: Map short reads for polishing
# ===========================================================================
log_step "L3: bwa short-read mapping for polishing"
if [ ! -f "${OUT}/polished/aligned.sam" ]; then
    bwa index "${OUT}/assembly/assembly.fasta"
    bwa mem -t ${THREADS} "${OUT}/assembly/assembly.fasta" \
            "${OUT}/trimmed/R1.fastq" "${OUT}/trimmed/R2.fastq" \
            > "${OUT}/polished/aligned.sam"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4 LEFT: Assembly QC
# ===========================================================================
log_step "L4-LEFT: QUAST"
if [ ! -f "${OUT}/qc/quast/report.tsv" ]; then
    quast "${OUT}/assembly/assembly.fasta" -r "${REFERENCE}" -o "${OUT}/qc/quast" -t ${THREADS}
else echo "Skipping (exists)"; fi

log_step "L4-LEFT: BUSCO"
if [ ! -d "${OUT}/qc/busco" ]; then
    busco -i "${OUT}/assembly/assembly.fasta" -o busco --out_path "${OUT}/qc" \
          -l bacteria_odb10 -m genome -c ${THREADS} --force
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4 RIGHT: Polypolish (short-read polishing of hybrid assembly)
# ===========================================================================
log_step "L4-RIGHT: Polypolish"
if [ ! -f "${OUT}/polished/polished.fasta" ]; then
    polypolish filter --in1 "${OUT}/polished/aligned.sam" --in2 "${OUT}/polished/aligned.sam" \
               --out1 "${OUT}/polished/filtered1.sam" --out2 "${OUT}/polished/filtered2.sam" 2>&1 || true
    polypolish polish "${OUT}/assembly/assembly.fasta" \
               "${OUT}/polished/aligned.sam" > "${OUT}/polished/polished.fasta" 2>&1 || {
        echo "WARNING: Polypolish failed, using unpolished assembly"
        cp "${OUT}/assembly/assembly.fasta" "${OUT}/polished/polished.fasta"
    }
else echo "Skipping (exists)"; fi
FINAL="${OUT}/polished/polished.fasta"

# ===========================================================================
# L5: Annotation + characterization
# ===========================================================================
log_step "L5: Prokka annotation"
if [ ! -f "${OUT}/prokka/HYBRID.gff" ]; then
    prokka "${FINAL}" --outdir "${OUT}/prokka" --prefix HYBRID \
           --cpus ${THREADS} --kingdom Bacteria --force
else echo "Skipping (exists)"; fi

log_step "L5: abricate AMR"
if [ ! -f "${OUT}/amr/abricate_card.tsv" ]; then
    abricate "${FINAL}" --db card --minid 80 --mincov 60 > "${OUT}/amr/abricate_card.tsv"
else echo "Skipping (exists)"; fi

log_step "L5: mlst typing"
if [ ! -f "${OUT}/mlst_out/mlst.tsv" ]; then
    mlst "${FINAL}" > "${OUT}/mlst_out/mlst.tsv"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L6: MERGE
# ===========================================================================
log_step "L6-MERGE: Building results"

TOTAL_LEN=$(grep "^Total length" "${OUT}/qc/quast/report.tsv" | head -1 | cut -f2)
NUM_CONTIGS=$(grep "^# contigs " "${OUT}/qc/quast/report.tsv" | head -1 | cut -f2)
N50=$(grep "^N50" "${OUT}/qc/quast/report.tsv" | cut -f2)
GC=$(grep "^GC" "${OUT}/qc/quast/report.tsv" | cut -f2)
LARGEST=$(grep "^Largest contig" "${OUT}/qc/quast/report.tsv" | cut -f2)
GENOME_FRAC=$(grep "^Genome fraction" "${OUT}/qc/quast/report.tsv" | cut -f2)

BUSCO_SUM=$(grep "C:" "${OUT}/qc/busco/short_summary.specific.bacteria_odb10.busco.txt" 2>/dev/null \
    | head -1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || echo "N/A")

CDS=$(grep "^CDS" "${OUT}/prokka/HYBRID.txt" | awk '{print $2}')
TRNA=$(grep "^tRNA" "${OUT}/prokka/HYBRID.txt" | awk '{print $2}')
RRNA=$(grep "^rRNA" "${OUT}/prokka/HYBRID.txt" | awk '{print $2}')

MLST_SCHEME=$(cut -f2 "${OUT}/mlst_out/mlst.tsv")
MLST_ST=$(cut -f3 "${OUT}/mlst_out/mlst.tsv")
AMR_COUNT=$(tail -n +2 "${OUT}/amr/abricate_card.tsv" 2>/dev/null | wc -l | tr -d ' ')

# Circularity from Unicycler
CIRCULAR=$(grep -c "circular=true" "${OUT}/assembly/assembly.fasta" || echo "0")

cat > "${RES}/hybrid_assembly_report.csv" << CSVEOF
metric,value
assembly_length,${TOTAL_LEN}
num_contigs,${NUM_CONTIGS}
n50,${N50}
gc_content,${GC}
largest_contig,${LARGEST}
genome_fraction,${GENOME_FRAC}
circular_contigs,${CIRCULAR}
completeness,${BUSCO_SUM}
mlst_scheme,${MLST_SCHEME}
mlst_sequence_type,${MLST_ST}
cds_count,${CDS}
trna_count,${TRNA}
rrna_count,${RRNA}
amr_genes,${AMR_COUNT}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/hybrid_assembly_report.csv"
echo ""
ls -lh "${RES}/"
