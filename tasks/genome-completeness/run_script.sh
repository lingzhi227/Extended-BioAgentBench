#!/bin/bash
set -e

# =============================================================================
# Task 29: Genome Completeness and Contamination Assessment
#
# DAG (depth 5):
# L0: assembled genome
# L1: ├── BUSCO (single-copy orthologs)
#     ├── QUAST (assembly metrics)
#     └── Prokka (gene prediction for stats)
# L2: ├── seqkit stats (basic stats)
#     └── parse BUSCO/QUAST/Prokka outputs
# L3: checkv (if viral) OR manual assessment
# L4: MERGE (comprehensive quality report)
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${DATA}/genome.fna"
REFERENCE="${REF}/reference.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{busco_out,quast,prokka,stats} "${RES}"

# L1-A: BUSCO
log_step "L1-A: BUSCO"
if [ ! -d "${OUT}/busco_out/busco" ]; then
    busco -i "${GENOME}" -o busco --out_path "${OUT}/busco_out" \
          -l bacteria_odb10 -m genome -c ${THREADS} --force
else echo "Skipping (exists)"; fi

# L1-B: QUAST with reference
log_step "L1-B: QUAST"
if [ ! -f "${OUT}/quast/report.tsv" ]; then
    quast "${GENOME}" -r "${REFERENCE}" -o "${OUT}/quast" -t ${THREADS}
else echo "Skipping (exists)"; fi

# L1-C: Prokka
log_step "L1-C: Prokka"
if [ ! -f "${OUT}/prokka/genome.gff" ]; then
    prokka "${GENOME}" --outdir "${OUT}/prokka" --prefix genome \
           --cpus ${THREADS} --kingdom Bacteria --force
else echo "Skipping (exists)"; fi

# L2: seqkit + parse
log_step "L2: Parse results"

seqkit stats "${GENOME}" -T > "${OUT}/stats/seqkit.tsv"

# BUSCO
BUSCO_SUM=$(grep "C:" "${OUT}/busco_out/busco/short_summary.specific.bacteria_odb10.busco.txt" 2>/dev/null \
    | head -1 | sed 's/^[[:space:]]*//' || echo "N/A")
BUSCO_COMPLETE=$(echo "$BUSCO_SUM" | grep -oP 'C:[\d.]+' | tr -d 'C:')
BUSCO_FRAGMENTED=$(echo "$BUSCO_SUM" | grep -oP 'F:[\d.]+' | tr -d 'F:')
BUSCO_MISSING=$(echo "$BUSCO_SUM" | grep -oP 'M:[\d.]+' | tr -d 'M:')

# QUAST
TOTAL_LEN=$(grep "^Total length" "${OUT}/quast/report.tsv" | head -1 | cut -f2)
NUM_CONTIGS=$(grep "^# contigs " "${OUT}/quast/report.tsv" | head -1 | cut -f2)
N50=$(grep "^N50" "${OUT}/quast/report.tsv" | cut -f2)
GC=$(grep "^GC" "${OUT}/quast/report.tsv" | cut -f2)
LARGEST=$(grep "^Largest contig" "${OUT}/quast/report.tsv" | cut -f2)
GF=$(grep "^Genome fraction" "${OUT}/quast/report.tsv" | cut -f2)
MISASSEMBLIES=$(grep "^# misassemblies" "${OUT}/quast/report.tsv" | cut -f2)
DUPLICATION=$(grep "^Duplication ratio" "${OUT}/quast/report.tsv" | cut -f2)

# Prokka
CDS=$(grep "^CDS" "${OUT}/prokka/genome.txt" | awk '{print $2}')
TRNA=$(grep "^tRNA" "${OUT}/prokka/genome.txt" | awk '{print $2}')
RRNA=$(grep "^rRNA" "${OUT}/prokka/genome.txt" | awk '{print $2}')
CODING_DENSITY=$(python3 -c "
cds_len = 0
for line in open('${OUT}/prokka/genome.gff'):
    if '\tCDS\t' in line:
        parts = line.split('\t')
        cds_len += int(parts[4]) - int(parts[3]) + 1
print(f'{cds_len/${TOTAL_LEN}*100:.1f}')
" 2>/dev/null || echo "N/A")

# L4: MERGE
log_step "L4-MERGE"

cat > "${RES}/completeness_report.csv" << CSVEOF
metric,value
total_length,${TOTAL_LEN}
num_contigs,${NUM_CONTIGS}
n50,${N50}
gc_content,${GC}
largest_contig,${LARGEST}
genome_fraction,${GF}
misassemblies,${MISASSEMBLIES}
duplication_ratio,${DUPLICATION}
busco_complete_pct,${BUSCO_COMPLETE}
busco_fragmented_pct,${BUSCO_FRAGMENTED}
busco_missing_pct,${BUSCO_MISSING}
cds_count,${CDS}
trna_count,${TRNA}
rrna_count,${RRNA}
coding_density_pct,${CODING_DENSITY}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/completeness_report.csv"
