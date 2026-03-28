#!/bin/bash
set -e

# =============================================================================
# Task 16: Long-read de novo Assembly and Polishing
#
# DAG (depth 7, linear chain with late fan-out):
#
# L0: raw nanopore reads
# L1: NanoPlot (read QC)
# L2: Filtlong (quality/length filter)
# L3: Flye (de novo assembly)
# L4: minimap2+samtools (map reads → sorted BAM for polishing)
# L5: Medaka (consensus polishing, round 1)
#     ├──────────────────────────────────────┐
# L6: QUAST (assembly QC vs ref)          Prokka (annotation)
#     │                                      │
#     ├── BUSCO (completeness)            abricate (AMR)
#     │                                      │
# L7: MERGE ─────────────────────────────────┘
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

READS="${DATA}/barcode10.fastq.gz"
REFERENCE="${REF}/salmonella_ref.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{nanoplot,filtered,assembly,mapping,polished,qc,prokka,amr} "${RES}"

# ===========================================================================
# L1: Read QC with NanoPlot
# ===========================================================================
log_step "L1: NanoPlot read QC"
if [ ! -f "${OUT}/nanoplot/NanoStats.txt" ]; then
    NanoPlot --fastq "${READS}" -o "${OUT}/nanoplot" -t ${THREADS} --plots dot
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2: Quality filtering with Filtlong
# ===========================================================================
log_step "L2: Filtlong quality filtering"
if [ ! -f "${OUT}/filtered/filtered.fastq.gz" ]; then
    filtlong --min_length 200 "${READS}" | gzip > "${OUT}/filtered/filtered.fastq.gz"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3: De novo assembly with Flye
# ===========================================================================
log_step "L3: Flye assembly"
if [ ! -f "${OUT}/assembly/assembly.fasta" ]; then
    flye --nano-raw "${OUT}/filtered/filtered.fastq.gz" \
         --out-dir "${OUT}/assembly" \
         --threads ${THREADS} --genome-size 5m
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4: Map reads back to assembly for polishing
# ===========================================================================
log_step "L4: minimap2 read mapping"
if [ ! -f "${OUT}/mapping/reads2assembly.bam" ]; then
    minimap2 -ax map-ont -t ${THREADS} \
             "${OUT}/assembly/assembly.fasta" "${OUT}/filtered/filtered.fastq.gz" \
        | samtools sort -@ ${THREADS} -o "${OUT}/mapping/reads2assembly.bam"
    samtools index "${OUT}/mapping/reads2assembly.bam"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L5: Consensus polishing with Medaka
# ===========================================================================
log_step "L5: Medaka polishing"
if [ ! -f "${OUT}/polished/consensus.fasta" ]; then
    medaka_consensus -i "${OUT}/filtered/filtered.fastq.gz" \
                     -d "${OUT}/assembly/assembly.fasta" \
                     -o "${OUT}/polished" \
                     -t ${THREADS} -m r1041_e82_400bps_sup_v5.0.0 2>&1 || {
        echo "WARNING: Medaka failed with specified model, trying default"
        medaka_consensus -i "${OUT}/filtered/filtered.fastq.gz" \
                         -d "${OUT}/assembly/assembly.fasta" \
                         -o "${OUT}/polished" \
                         -t ${THREADS} 2>&1 || {
            echo "WARNING: Medaka failed, using unpolished assembly"
            cp "${OUT}/assembly/assembly.fasta" "${OUT}/polished/consensus.fasta"
        }
    }
else echo "Skipping (exists)"; fi
FINAL="${OUT}/polished/consensus.fasta"

# ===========================================================================
# L6 LEFT: Assembly QC
# ===========================================================================
log_step "L6-LEFT: QUAST assembly QC"
if [ ! -f "${OUT}/qc/quast/report.tsv" ]; then
    quast "${FINAL}" -r "${REFERENCE}" -o "${OUT}/qc/quast" -t ${THREADS}
else echo "Skipping (exists)"; fi

log_step "L6-LEFT: BUSCO completeness"
if [ ! -d "${OUT}/qc/busco" ]; then
    busco -i "${FINAL}" -o busco --out_path "${OUT}/qc" \
          -l bacteria_odb10 -m genome -c ${THREADS} --force
else echo "Skipping (exists)"; fi

# ===========================================================================
# L6 RIGHT: Annotation + AMR
# ===========================================================================
log_step "L6-RIGHT: Prokka annotation"
if [ ! -f "${OUT}/prokka/SALM.gff" ]; then
    prokka "${FINAL}" --outdir "${OUT}/prokka" --prefix SALM \
           --cpus ${THREADS} --kingdom Bacteria --genus Salmonella --force
else echo "Skipping (exists)"; fi

log_step "L6-RIGHT: abricate AMR detection"
if [ ! -f "${OUT}/amr/abricate_card.tsv" ]; then
    abricate "${FINAL}" --db card --minid 80 --mincov 60 > "${OUT}/amr/abricate_card.tsv"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L7: MERGE
# ===========================================================================
log_step "L7-MERGE: Building results"

# Read QC stats
MEAN_LEN=$(grep "Mean read length" "${OUT}/nanoplot/NanoStats.txt" | awk '{print $NF}' | tr -d ',')
MEAN_QUAL=$(grep "Mean read quality" "${OUT}/nanoplot/NanoStats.txt" | awk '{print $NF}')
TOTAL_BASES=$(grep "Total bases" "${OUT}/nanoplot/NanoStats.txt" | awk '{print $NF}' | tr -d ',')
NUM_READS=$(grep "Number of reads" "${OUT}/nanoplot/NanoStats.txt" | head -1 | awk '{print $NF}' | tr -d ',')

# Assembly stats from QUAST
TOTAL_LEN=$(grep "^Total length" "${OUT}/qc/quast/report.tsv" | head -1 | cut -f2)
NUM_CONTIGS=$(grep "^# contigs " "${OUT}/qc/quast/report.tsv" | head -1 | cut -f2)
N50=$(grep "^N50" "${OUT}/qc/quast/report.tsv" | cut -f2)
GC=$(grep "^GC" "${OUT}/qc/quast/report.tsv" | cut -f2)
LARGEST=$(grep "^Largest contig" "${OUT}/qc/quast/report.tsv" | cut -f2)
GENOME_FRAC=$(grep "^Genome fraction" "${OUT}/qc/quast/report.tsv" | cut -f2)
MISASSEMBLIES=$(grep "^# misassemblies" "${OUT}/qc/quast/report.tsv" | cut -f2)

# BUSCO
BUSCO_SUM=$(grep "C:" "${OUT}/qc/busco/short_summary.specific.bacteria_odb10.busco.txt" 2>/dev/null \
    | head -1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || echo "N/A")

# Prokka
CDS=$(grep "^CDS" "${OUT}/prokka/SALM.txt" | awk '{print $2}')
TRNA=$(grep "^tRNA" "${OUT}/prokka/SALM.txt" | awk '{print $2}')
RRNA=$(grep "^rRNA" "${OUT}/prokka/SALM.txt" | awk '{print $2}')

# AMR
AMR_COUNT=$(tail -n +2 "${OUT}/amr/abricate_card.tsv" 2>/dev/null | wc -l | tr -d ' ')

cat > "${RES}/longread_assembly_report.csv" << CSVEOF
metric,value
num_reads,${NUM_READS}
mean_read_length,${MEAN_LEN}
mean_read_quality,${MEAN_QUAL}
total_bases,${TOTAL_BASES}
assembly_length,${TOTAL_LEN}
num_contigs,${NUM_CONTIGS}
n50,${N50}
gc_content,${GC}
largest_contig,${LARGEST}
genome_fraction,${GENOME_FRAC}
misassemblies,${MISASSEMBLIES}
completeness,${BUSCO_SUM}
cds_count,${CDS}
trna_count,${TRNA}
rrna_count,${RRNA}
amr_genes,${AMR_COUNT}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/longread_assembly_report.csv"
echo ""
ls -lh "${RES}/"
