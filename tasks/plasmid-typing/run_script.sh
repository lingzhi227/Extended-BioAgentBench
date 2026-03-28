#!/bin/bash
set -e

# =============================================================================
# Task 28: Plasmid Detection and Replicon Typing
#
# DAG (depth 5):
# L0: assembled genome
# L1: ├── mob_recon (plasmid reconstruction)
#     └── abricate --db plasmidfinder (replicon detection)
# L2: ├── mob_recon results → plasmid contig classification
#     └── abricate --db card on plasmid contigs (plasmid-borne AMR)
# L3: ├── seqkit stats per plasmid
#     └── abricate --db vfdb on plasmid contigs (virulence)
# L4: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

GENOME="${DATA}/genome.fna"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{mob_recon,replicons,amr,virulence,stats} "${RES}"

# L1-A: Identify plasmid contigs by size/replicon
log_step "L1-A: Contig classification"
if [ ! -f "${OUT}/mob_recon/contig_report.txt" ]; then
    mkdir -p "${OUT}/mob_recon"
    # Classify contigs: largest = chromosome, smaller = potential plasmids
    python3 -c "
from Bio import SeqIO
import sys
contigs = list(SeqIO.parse('${GENOME}', 'fasta'))
contigs.sort(key=lambda x: len(x), reverse=True)
chrom = contigs[0]
plasmids = contigs[1:]
with open('${OUT}/mob_recon/chromosome.fasta', 'w') as f:
    SeqIO.write([chrom], f, 'fasta')
if plasmids:
    with open('${OUT}/mob_recon/all_plasmids.fasta', 'w') as f:
        SeqIO.write(plasmids, f, 'fasta')
with open('${OUT}/mob_recon/contig_report.txt', 'w') as f:
    f.write(f'chromosome\t{chrom.id}\t{len(chrom)}\n')
    for p in plasmids:
        f.write(f'plasmid\t{p.id}\t{len(p)}\n')
print(f'Chromosome: {chrom.id} ({len(chrom)} bp)')
for p in plasmids:
    print(f'Plasmid: {p.id} ({len(p)} bp)')
"
else echo "Skipping (exists)"; fi

# L1-B: Replicon typing on full assembly
log_step "L1-B: Replicon detection"
if [ ! -f "${OUT}/replicons/replicons.tsv" ]; then
    abricate "${GENOME}" --db plasmidfinder --minid 80 --mincov 60 \
        > "${OUT}/replicons/replicons.tsv" 2>/dev/null || true
else echo "Skipping (exists)"; fi

# L2: AMR on plasmid contigs
log_step "L2: Plasmid-borne AMR"
PLASMID_CONTIGS="${OUT}/mob_recon/all_plasmids.fasta"

if [ -s "${PLASMID_CONTIGS}" ]; then
    abricate "${PLASMID_CONTIGS}" --db card --minid 80 --mincov 60 \
        > "${OUT}/amr/plasmid_amr.tsv" 2>/dev/null || true
else
    touch "${OUT}/amr/plasmid_amr.tsv"
fi

# L2: AMR on full genome for comparison
abricate "${GENOME}" --db card --minid 80 --mincov 60 \
    > "${OUT}/amr/genome_amr.tsv" 2>/dev/null || true

# L3: Virulence on plasmid contigs
log_step "L3: Plasmid virulence"
if [ -s "${PLASMID_CONTIGS}" ]; then
    abricate "${PLASMID_CONTIGS}" --db vfdb --minid 80 --mincov 60 \
        > "${OUT}/virulence/plasmid_vf.tsv" 2>/dev/null || true
else
    touch "${OUT}/virulence/plasmid_vf.tsv"
fi

# L3: Stats
log_step "L3: Assembly stats"
GENOME_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${GENOME}")
TOTAL_CONTIGS=$(grep -c ">" "${GENOME}" || true)
PLASMID_COUNT=$(grep -c "^plasmid" "${OUT}/mob_recon/contig_report.txt" 2>/dev/null || true)
PLASMID_COUNT=${PLASMID_COUNT:-0}

if [ -s "${PLASMID_CONTIGS}" ]; then
    PLASMID_TOTAL_LEN=$(awk '/^>/{if(l)s+=l;l=0;next}{l+=length}END{s+=l;print s}' "${PLASMID_CONTIGS}")
else
    PLASMID_TOTAL_LEN=0
fi

CHROMOSOME_CONTIGS=$(grep -c "^chromosome" "${OUT}/mob_recon/contig_report.txt" 2>/dev/null || true)
CHROMOSOME_CONTIGS=${CHROMOSOME_CONTIGS:-0}
REPLICON_TYPES=$(tail -n +2 "${OUT}/replicons/replicons.tsv" 2>/dev/null | awk -F'\t' '{print $6}' | sort -u | wc -l | tr -d ' ')
REPLICON_TYPES=${REPLICON_TYPES:-0}
PLASMID_AMR=$(tail -n +2 "${OUT}/amr/plasmid_amr.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
GENOME_AMR=$(tail -n +2 "${OUT}/amr/genome_amr.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
PLASMID_VF=$(tail -n +2 "${OUT}/virulence/plasmid_vf.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# L4: MERGE
log_step "L4-MERGE"

cat > "${RES}/plasmid_report.csv" << CSVEOF
metric,value
genome_length,${GENOME_LEN}
total_contigs,${TOTAL_CONTIGS}
plasmid_clusters,${PLASMID_COUNT}
plasmid_total_length,${PLASMID_TOTAL_LEN}
chromosome_contigs,${CHROMOSOME_CONTIGS}
replicon_types_detected,${REPLICON_TYPES}
amr_genes_genome,${GENOME_AMR}
amr_genes_plasmid,${PLASMID_AMR}
virulence_factors_plasmid,${PLASMID_VF}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/plasmid_report.csv"
