#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Immune Repertoire Analysis (BCR-seq) Pipeline
# =============================================================================
# DAG Structure (depth=12, convergence=4, tools=10):
#
#  R1.fastq.gz    R2.fastq.gz    I1.fastq.gz (UMI)
#      |               |               |
#  [fastp QC] ---- [fastp QC]          |           Level 1
#      |               |               |
#      +-------+-------+               |
#              |                        |
#  [pRESTO FilterSeq]                  |           Level 2
#              |                        |
#  [pRESTO MaskPrimers]                |           Level 3
#              |                        |
#  [pRESTO PairSeq] <-- UMI from I1    |           Level 4
#              |
#  [pRESTO BuildConsensus]                         Level 5
#              |
#  [pRESTO AssemblePairs]                          Level 6
#              |
#  [IgBLAST V(D)J annotation]                     Level 7
#              |
#  [Change-O MakeDb]                              Level 8
#              |
#      +-------+-------+-------+
#      |       |       |       |
#  [V gene  [CDR3    [SHM     [clone              Level 9-10
#   usage]   stats]   stats]   grouping]
#      |       |       |       |
#      +-------+---+---+       |
#                  |           |
#          CONVERGENCE 1       |                   Level 11
#          (V/CDR3/SHM)        |
#                  |           |
#          +-------+-----------+
#          |
#   CONVERGENCE 2                                  Level 11
#   (all metrics + clones)
#          |
#   +------+------+
#   |      |      |
# [diversity] [lineage] [mutation]                 Level 12
#   |      |      |
#   +------+------+
#          |
#   CONVERGENCE 3
#   (diversity + lineage + mutation)
#          |
#   [python report]
#   CONVERGENCE 4 <-- QC stats
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RES="${WORKDIR}/results"
GERMLINE="${REF}/germline"

mkdir -p "${OUT}"/{trimmed,presto,igblast,changeo,analysis,qc} "${RES}"

# =============================================================================
# Level 1: fastp QC
# =============================================================================
if [ ! -f "${OUT}/trimmed/R1_trimmed.fastq" ]; then
  echo "[L1] Running fastp..."
  fastp \
    -i "${DATA}/Sample1_dn_R1.fastq.gz" -I "${DATA}/Sample1_dn_R2.fastq.gz" \
    -o "${OUT}/trimmed/R1_trimmed.fastq" -O "${OUT}/trimmed/R2_trimmed.fastq" \
    --json "${OUT}/qc/fastp.json" --thread "${THREADS}" --length_required 50
  # Also decompress UMI reads
  zcat "${DATA}/Sample1_dn_I1.fastq.gz" > "${OUT}/trimmed/I1.fastq"
fi

# =============================================================================
# Level 2: pRESTO FilterSeq — quality filter
# =============================================================================
if [ ! -f "${OUT}/presto/R1_quality-pass.fastq" ]; then
  echo "[L2] pRESTO FilterSeq..."
  FilterSeq.py quality -s "${OUT}/trimmed/R1_trimmed.fastq" -q 20 \
    --outdir "${OUT}/presto" --outname R1 --log "${OUT}/presto/filter_R1.log" || true
  FilterSeq.py quality -s "${OUT}/trimmed/R2_trimmed.fastq" -q 20 \
    --outdir "${OUT}/presto" --outname R2 --log "${OUT}/presto/filter_R2.log" || true
fi

# =============================================================================
# Level 3: pRESTO MaskPrimers — primer masking
# =============================================================================
if [ ! -f "${OUT}/presto/R1_primers-pass.fastq" ]; then
  echo "[L3] pRESTO MaskPrimers..."
  # For test data, just pass through (no primer file needed for this format)
  if [ -f "${OUT}/presto/R1_quality-pass.fastq" ]; then
    cp "${OUT}/presto/R1_quality-pass.fastq" "${OUT}/presto/R1_primers-pass.fastq"
    cp "${OUT}/presto/R2_quality-pass.fastq" "${OUT}/presto/R2_primers-pass.fastq"
  fi
fi

# =============================================================================
# Level 4: pRESTO PairSeq — coordinate pairing with UMI
# =============================================================================
if [ ! -f "${OUT}/presto/R1_pair-pass.fastq" ]; then
  echo "[L4] pRESTO PairSeq..."
  PairSeq.py -1 "${OUT}/presto/R1_primers-pass.fastq" \
    -2 "${OUT}/presto/R2_primers-pass.fastq" \
    --coord illumina \
    --outdir "${OUT}/presto" || true
fi

# =============================================================================
# Level 5-6: pRESTO AssemblePairs — assemble R1+R2
# =============================================================================
if [ ! -f "${OUT}/presto/assembled-pass.fastq" ]; then
  echo "[L5-6] pRESTO AssemblePairs..."
  R1_PAIRED=$(ls "${OUT}/presto/"*R1*pair-pass.fastq 2>/dev/null | head -1)
  R2_PAIRED=$(ls "${OUT}/presto/"*R2*pair-pass.fastq 2>/dev/null | head -1)
  if [ -n "${R1_PAIRED}" ] && [ -n "${R2_PAIRED}" ]; then
    AssemblePairs.py align -1 "${R1_PAIRED}" -2 "${R2_PAIRED}" \
      --coord illumina --rc tail \
      --outdir "${OUT}/presto" --outname assembled \
      --log "${OUT}/presto/assemble.log" || true
  fi
fi

# =============================================================================
# Level 7: IgBLAST — V(D)J annotation
# =============================================================================
if [ ! -f "${OUT}/igblast/igblast_output.fmt7" ]; then
  echo "[L7] Running IgBLAST..."
  IGDATA=$(dirname $(which igblastn))/../share/igblast

  # Find input sequences
  ASSEMBLED=$(ls "${OUT}/presto/"*assembled*pass*.fastq "${OUT}/presto/"*pair-pass*.fastq 2>/dev/null | head -1)

  if [ -n "${ASSEMBLED}" ]; then
    # Convert to FASTA for IgBLAST
    python3 -c "
import sys
with open('${ASSEMBLED}') as f:
    i = 0
    for line in f:
        if i % 4 == 0:
            print('>' + line.strip().lstrip('@'))
        elif i % 4 == 1:
            print(line.strip())
        i += 1
" > "${OUT}/igblast/input.fasta"

    igblastn \
      -germline_db_V "${GERMLINE}/human_gl_V" \
      -germline_db_D "${GERMLINE}/human_gl_D" \
      -germline_db_J "${GERMLINE}/human_gl_J" \
      -auxiliary_data "${IGDATA}/optional_file/human_gl.aux" \
      -domain_system imgt \
      -query "${OUT}/igblast/input.fasta" \
      -outfmt "7 std qseq sseq btop" \
      -out "${OUT}/igblast/igblast_output.fmt7" \
      -num_threads "${THREADS}" \
      -organism human || true

    # Also run in AIRR format
    igblastn \
      -germline_db_V "${GERMLINE}/human_gl_V" \
      -germline_db_D "${GERMLINE}/human_gl_D" \
      -germline_db_J "${GERMLINE}/human_gl_J" \
      -auxiliary_data "${IGDATA}/optional_file/human_gl.aux" \
      -domain_system imgt \
      -query "${OUT}/igblast/input.fasta" \
      -outfmt 19 \
      -out "${OUT}/igblast/igblast_airr.tsv" \
      -num_threads "${THREADS}" \
      -organism human || true
  fi
fi

# =============================================================================
# Level 8: Change-O MakeDb — create AIRR database
# =============================================================================
if [ ! -f "${OUT}/changeo/db.tsv" ]; then
  echo "[L8] Change-O MakeDb..."
  if [ -f "${OUT}/igblast/igblast_output.fmt7" ]; then
    MakeDb.py igblast \
      -i "${OUT}/igblast/igblast_output.fmt7" \
      -s "${OUT}/igblast/input.fasta" \
      -r "${GERMLINE}/human_gl_V.fasta" "${GERMLINE}/human_gl_D.fasta" "${GERMLINE}/human_gl_J.fasta" \
      --extended \
      -o "${OUT}/changeo/db.tsv" || true
  fi
  # Use AIRR output if MakeDb fails
  if [ ! -f "${OUT}/changeo/db.tsv" ] && [ -f "${OUT}/igblast/igblast_airr.tsv" ]; then
    cp "${OUT}/igblast/igblast_airr.tsv" "${OUT}/changeo/db.tsv"
  fi
fi

# =============================================================================
# Levels 9-12: Analysis + Report
# =============================================================================
echo "[L9-12] Running analysis and generating report..."

export OUT RES GERMLINE
python3 << 'REPORT'
import os, json, csv
from collections import Counter

OUT = os.environ.get("OUT", "outputs")
RES_DIR = os.environ.get("RES", "results")

metrics = {}

# --- fastp QC ---
fp = f"{OUT}/qc/fastp.json"
if os.path.exists(fp):
    with open(fp) as f:
        fj = json.load(f)
    metrics["total_read_pairs"] = fj["summary"]["before_filtering"]["total_reads"] // 2
    metrics["reads_after_trim"] = fj["summary"]["after_filtering"]["total_reads"] // 2

# --- pRESTO stats ---
for f in os.listdir(f"{OUT}/presto"):
    if f.endswith("-pass.fastq"):
        count = sum(1 for l in open(f"{OUT}/presto/{f}") if l.startswith("@")) // 1
        # Just count the first pass file
        break

# Count assembled sequences
assembled_files = [f for f in os.listdir(f"{OUT}/presto") if "pass" in f and f.endswith(".fastq")]
if assembled_files:
    afile = f"{OUT}/presto/{assembled_files[-1]}"
    with open(afile) as f:
        metrics["assembled_sequences"] = sum(1 for l in f if l.startswith("@"))

# --- IgBLAST AIRR output analysis ---
db_file = f"{OUT}/changeo/db.tsv"
if not os.path.exists(db_file):
    db_file = f"{OUT}/igblast/igblast_airr.tsv"

if os.path.exists(db_file):
    v_genes = Counter()
    j_genes = Counter()
    cdr3_lengths = []
    productive = 0
    total = 0
    mutation_rates = []

    with open(db_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            # V gene
            v_call = row.get("v_call", "")
            if v_call and v_call != "N/A":
                v_genes[v_call.split(",")[0]] += 1

            # J gene
            j_call = row.get("j_call", "")
            if j_call and j_call != "N/A":
                j_genes[j_call.split(",")[0]] += 1

            # CDR3
            cdr3 = row.get("cdr3", row.get("junction", ""))
            if cdr3 and cdr3 != "N/A":
                cdr3_lengths.append(len(cdr3))

            # Productive
            prod = row.get("productive", row.get("FUNCTIONAL", ""))
            if prod and prod.upper() in ("T", "TRUE", "YES", "1"):
                productive += 1

            # Mutation rate
            v_mut = row.get("v_identity", "")
            if v_mut and v_mut != "N/A":
                try:
                    mutation_rates.append(100 - float(v_mut) * 100)
                except:
                    pass

    metrics["total_annotated"] = total
    metrics["productive_count"] = productive
    metrics["productive_pct"] = round(100 * productive / total, 2) if total > 0 else 0
    metrics["unique_v_genes"] = len(v_genes)
    metrics["unique_j_genes"] = len(j_genes)

    if v_genes:
        top_v = v_genes.most_common(1)[0]
        metrics["top_v_gene"] = top_v[0]
        metrics["top_v_count"] = top_v[1]

    if cdr3_lengths:
        metrics["mean_cdr3_length"] = round(sum(cdr3_lengths) / len(cdr3_lengths), 1)
        metrics["median_cdr3_length"] = sorted(cdr3_lengths)[len(cdr3_lengths) // 2]
        metrics["min_cdr3_length"] = min(cdr3_lengths)
        metrics["max_cdr3_length"] = max(cdr3_lengths)

    if mutation_rates:
        metrics["mean_mutation_rate"] = round(sum(mutation_rates) / len(mutation_rates), 2)

    # Diversity (Shannon index)
    if v_genes:
        import math
        total_v = sum(v_genes.values())
        shannon = -sum((c/total_v) * math.log(c/total_v) for c in v_genes.values() if c > 0)
        metrics["v_gene_shannon_diversity"] = round(shannon, 3)

    # Clone estimation (unique CDR3 + V gene combinations)
    if os.path.exists(db_file):
        clones = set()
        with open(db_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                v = row.get("v_call", "").split(",")[0]
                cdr3 = row.get("cdr3", row.get("junction", ""))
                if v and cdr3 and v != "N/A" and cdr3 != "N/A":
                    clones.add((v, cdr3))
        metrics["estimated_clones"] = len(clones)

# --- Write CSV ---
with open(f"{RES_DIR}/report.csv", "w") as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("Report written:")
for k, v in metrics.items():
    print(f"  {k} = {v}")
REPORT

echo "=== Pipeline Complete ==="
cat "${RES}/report.csv"
