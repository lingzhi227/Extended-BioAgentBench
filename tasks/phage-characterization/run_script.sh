#!/bin/bash
set -e

# =============================================================================
# Task 21: Phage Genome Characterization
#
# DAG (depth 5, fan-out from assembly):
#
# L0: reads
# L1: fastp (trim)
# L2: shovill (assemble)
#     ├────────────────────────────────────────┐
# L3: pharokka (phage annotation)       checkv (quality)
#     │                                        │
# L4: functional summary                 completeness report
#     │                                        │
#     ├── abricate VFDB (virulence)            │
#     └────────────────────────────────────────┘
# L5: MERGE
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${SCRIPT_DIR}/data"
REF="${SCRIPT_DIR}/reference"
OUT="${SCRIPT_DIR}/outputs"
RES="${SCRIPT_DIR}/results"

PHAROKKA_DB="/pscratch/sd/l/lingzhi/pharokka_db"
CHECKV_DB="/pscratch/sd/l/lingzhi/checkv_db/checkv-db-v1.5"

log_step() {
    echo "=================================================================="
    echo "STEP: $1"
    echo "$(date)"
    echo "=================================================================="
}

mkdir -p "${OUT}"/{trimmed,assembly,pharokka,checkv,vfdb} "${RES}"

# ===========================================================================
# L1: Trimming
# ===========================================================================
log_step "L1: fastp"
if [ ! -f "${OUT}/trimmed/R1.fastq.gz" ]; then
    fastp --in1 "${DATA}/reads_R1.fastq.gz" --in2 "${DATA}/reads_R2.fastq.gz" \
          --out1 "${OUT}/trimmed/R1.fastq.gz" --out2 "${OUT}/trimmed/R2.fastq.gz" \
          --detect_adapter_for_pe --thread ${THREADS} \
          --json "${OUT}/trimmed/fastp.json"
else echo "Skipping (exists)"; fi

# ===========================================================================
# L2: Assembly
# ===========================================================================
log_step "L2: MEGAHIT assembly"
if [ ! -f "${OUT}/assembly/final.contigs.fa" ]; then
    megahit -1 "${OUT}/trimmed/R1.fastq.gz" -2 "${OUT}/trimmed/R2.fastq.gz" \
            -o "${OUT}/assembly" -t ${THREADS} --min-contig-len 500
else echo "Skipping (exists)"; fi
CONTIGS="${OUT}/assembly/final.contigs.fa"

# ===========================================================================
# L3-A: Pharokka phage-specific annotation
# ===========================================================================
log_step "L3-A: pharokka annotation"
if [ ! -f "${OUT}/pharokka/pharokka_final.gff" ]; then
    pharokka.py -i "${CONTIGS}" -o "${OUT}/pharokka" \
                -d "${PHAROKKA_DB}" -t ${THREADS} --force 2>&1 || {
        echo "WARNING: pharokka failed, falling back to prokka --kingdom Viruses"
        prokka "${CONTIGS}" --outdir "${OUT}/pharokka" --prefix pharokka_final \
               --cpus ${THREADS} --kingdom Viruses --force 2>&1 || true
    }
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3-B: CheckV quality assessment
# ===========================================================================
log_step "L3-B: checkv quality"
if [ ! -f "${OUT}/checkv/quality_summary.tsv" ]; then
    if [ -d "${CHECKV_DB}" ]; then
        checkv end_to_end "${CONTIGS}" "${OUT}/checkv" -d "${CHECKV_DB}" -t ${THREADS} 2>&1 || true
    else
        echo "WARNING: CheckV database not found"
        mkdir -p "${OUT}/checkv"
        touch "${OUT}/checkv/quality_summary.tsv"
    fi
else echo "Skipping (exists)"; fi

# ===========================================================================
# L3-C: Virulence factor detection
# ===========================================================================
log_step "L3-C: abricate VFDB"
if [ ! -f "${OUT}/vfdb/vfdb_results.tsv" ]; then
    abricate "${CONTIGS}" --db vfdb --minid 80 --mincov 60 > "${OUT}/vfdb/vfdb_results.tsv" 2>/dev/null || true
else echo "Skipping (exists)"; fi

# ===========================================================================
# L4: Parse results
# ===========================================================================
log_step "L4: Parse annotation results"

# Pharokka stats
PHAROKKA_GFF="${OUT}/pharokka/pharokka.gff"
if [ -f "${PHAROKKA_GFF}" ]; then
    CDS=$(grep -c "CDS" "${PHAROKKA_GFF}" 2>/dev/null || echo "0")
    TRNA=$(grep -c "tRNA" "${OUT}/pharokka/pharokka_aragorn.gff" 2>/dev/null || true)
    TRNA=${TRNA:-0}
    # Count functional categories from pharokka GFF product annotations
    LYSIS=$(grep -ic "lysis\|holin\|endolysin\|spanin\|lysozyme" "${PHAROKKA_GFF}" 2>/dev/null || echo "0")
    LYSOGENY=$(grep -ic "integrase\|repressor\|lysogeny\|excisionase" "${PHAROKKA_GFF}" 2>/dev/null || echo "0")
    REPLICATION=$(grep -ic "helicase\|primase\|polymerase\|replication\|SSB" "${PHAROKKA_GFF}" 2>/dev/null || echo "0")
    STRUCTURAL=$(grep -ic "capsid\|tail\|baseplate\|head\|portal\|sheath\|fiber" "${PHAROKKA_GFF}" 2>/dev/null || echo "0")
    # Get top MASH hit for taxonomy
    TAXONOMY=$(tail -1 "${OUT}/pharokka/pharokka_top_hits_mash_inphared.tsv" 2>/dev/null | cut -f6 || echo "Unknown")
else
    CDS=0; TRNA=0; LYSIS=0; LYSOGENY=0; REPLICATION=0; STRUCTURAL=0; TAXONOMY="Unknown"
fi

# CheckV stats
if [ -s "${OUT}/checkv/quality_summary.tsv" ]; then
    QUALITY=$(tail -1 "${OUT}/checkv/quality_summary.tsv" | cut -f8 || echo "N/A")
    COMPLETENESS=$(tail -1 "${OUT}/checkv/quality_summary.tsv" | cut -f10 || echo "N/A")
    CONTAMINATION=$(tail -1 "${OUT}/checkv/quality_summary.tsv" | cut -f12 || echo "N/A")
else
    QUALITY="N/A"; COMPLETENESS="N/A"; CONTAMINATION="N/A"
fi

# Assembly stats
NUM_CONTIGS=$(grep -c ">" "${CONTIGS}" || echo "0")
TOTAL_LEN=$(awk '/^>/{if(l)print l; l=0; next}{l+=length}END{print l}' "${CONTIGS}" || echo "0")
GC=$(awk '/^>/{next}{s+=length; gc+=gsub(/[GCgc]/,"&")}END{printf "%.2f",gc/s*100}' "${CONTIGS}" || echo "0")

# Virulence
VF_COUNT=$(tail -n +2 "${OUT}/vfdb/vfdb_results.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# ===========================================================================
# L5: MERGE
# ===========================================================================
log_step "L5-MERGE"

cat > "${RES}/phage_report.csv" << CSVEOF
metric,value
genome_length,${TOTAL_LEN}
num_contigs,${NUM_CONTIGS}
gc_content,${GC}
closest_hit,${TAXONOMY}
cds_count,${CDS}
trna_count,${TRNA}
lysis_genes,${LYSIS}
lysogeny_genes,${LYSOGENY}
replication_genes,${REPLICATION}
structural_genes,${STRUCTURAL}
checkv_quality,${QUALITY}
checkv_completeness,${COMPLETENESS}
virulence_factors,${VF_COUNT}
CSVEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RES}/phage_report.csv"
echo ""
ls -lh "${RES}/"
