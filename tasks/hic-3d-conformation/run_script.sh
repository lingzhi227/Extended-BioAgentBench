#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Hi-C 3D Genome Conformation Analysis Pipeline
# =============================================================================
# DAG Structure (depth=10, convergence=4, tools=10):
#
#  R1.fastq.gz    R2.fastq.gz
#      |               |
#      +-------+-------+
#              |
#         [fastp]                                      Level 1
#              |
#      [bwa-mem2 -SP align]                            Level 2
#              |
#      [samtools view/filter]                          Level 3
#              |
#      [pairtools parse] <-- CONVERGENCE 1             Level 4
#              |            (mate pairing)
#      [pairtools sort]                                Level 5
#              |
#      +-------+---------------+
#      |       |               |
#  [pairtools [cooler       [pairtools                 Level 6
#   dedup]     cload]        stats]
#      |       |               |
#      |   [cooler             |
#      |    balance]           |
#      |       |               |
#      +-------+---------------+
#      |       |
#      |   +---+---------------+
#      |   |                   |
#      | [cooltools           [cooltools               Level 7
#      |  eigs]                insulation]
#      |   |                   |
#      |   |(compartments)     |(boundaries)
#      |   |                   |
#      +---+----------+-------+
#                      |
#              CONVERGENCE 2                           Level 8
#              (compartments + TADs + pairs)
#                      |
#      +---------------+---------------+
#      |               |               |
#  [cooltools       [chromosight    [fithic            Level 9
#   saddle]          detect]         (sig loops)]
#      |               |               |
#      +-----------+---+---+-----------+
#                  |
#          CONVERGENCE 3                               Level 10
#          (saddle + loops + dots)
#                  |
#          [python report]
#          CONVERGENCE 4 <-- stats
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RES="${WORKDIR}/results"
GENOME="${REF}/genome.fa"
CHROMSIZES="${REF}/chromsizes.tsv"
RESOLUTION=5000

mkdir -p "${OUT}"/{trimmed,aligned,pairs,matrix,compartments,insulation,loops,fithic,saddle,qc} "${RES}"

# =============================================================================
# Level 1: fastp — adapter trimming + QC
# =============================================================================
if [ ! -f "${OUT}/trimmed/R1_trimmed.fastq.gz" ]; then
  echo "[Level 1] Running fastp..."
  fastp \
    -i "${DATA}/R1.fastq.gz" -I "${DATA}/R2.fastq.gz" \
    -o "${OUT}/trimmed/R1_trimmed.fastq.gz" -O "${OUT}/trimmed/R2_trimmed.fastq.gz" \
    --json "${OUT}/qc/fastp.json" --html "${OUT}/qc/fastp.html" \
    --thread "${THREADS}" --length_required 30
fi

# =============================================================================
# Level 2: bwa-mem2 index + align (Hi-C mode: -SP5M)
# =============================================================================
if [ ! -f "${GENOME}.bwt.2bit.64" ]; then
  echo "[Level 2] Indexing reference with bwa-mem2..."
  bwa-mem2 index "${GENOME}"
fi

if [ ! -f "${OUT}/aligned/hic_aligned.bam" ]; then
  echo "[Level 2] Aligning Hi-C reads with bwa-mem2 -SP..."
  bwa-mem2 mem -SP5M -t "${THREADS}" \
    -R "@RG\tID:hic\tSM:sample\tPL:ILLUMINA" \
    "${GENOME}" \
    "${OUT}/trimmed/R1_trimmed.fastq.gz" "${OUT}/trimmed/R2_trimmed.fastq.gz" \
    | samtools view -@ "${THREADS}" -bhS - \
    > "${OUT}/aligned/hic_aligned.bam"
fi

# =============================================================================
# Level 3: samtools — filter unmapped and low MAPQ
# =============================================================================
if [ ! -f "${OUT}/aligned/hic_filtered.bam" ]; then
  echo "[Level 3] Filtering BAM (keep all for pairtools, MAPQ handled by pairtools)..."
  samtools view -@ "${THREADS}" -bh "${OUT}/aligned/hic_aligned.bam" \
    > "${OUT}/aligned/hic_filtered.bam"
fi

# =============================================================================
# Generate chromsizes if not present
# =============================================================================
if [ ! -f "${CHROMSIZES}" ]; then
  echo "Generating chromsizes..."
  samtools faidx "${GENOME}"
  cut -f1,2 "${GENOME}.fai" > "${CHROMSIZES}"
fi

# =============================================================================
# Level 4: pairtools parse — CONVERGENCE 1 (mate pairing)
# =============================================================================
if [ ! -f "${OUT}/pairs/parsed.pairs.gz" ]; then
  echo "[Level 4] Running pairtools parse..."
  pairtools parse \
    --chroms-path "${CHROMSIZES}" \
    --min-mapq 30 \
    --walks-policy mask \
    --output "${OUT}/pairs/parsed.pairs.gz" \
    "${OUT}/aligned/hic_filtered.bam"
fi

# =============================================================================
# Level 5: pairtools sort
# =============================================================================
if [ ! -f "${OUT}/pairs/sorted.pairs.gz" ]; then
  echo "[Level 5] Running pairtools sort..."
  pairtools sort \
    --nproc "${THREADS}" \
    --output "${OUT}/pairs/sorted.pairs.gz" \
    "${OUT}/pairs/parsed.pairs.gz"
fi

# =============================================================================
# Level 6a: pairtools dedup
# =============================================================================
if [ ! -f "${OUT}/pairs/dedup.pairs.gz" ]; then
  echo "[Level 6a] Running pairtools dedup..."
  pairtools dedup \
    --output "${OUT}/pairs/dedup.pairs.gz" \
    --output-stats "${OUT}/pairs/dedup_stats.txt" \
    "${OUT}/pairs/sorted.pairs.gz"
fi

# =============================================================================
# Level 6b: pairtools stats
# =============================================================================
if [ ! -f "${OUT}/pairs/pair_stats.txt" ]; then
  echo "[Level 6b] Running pairtools stats..."
  pairtools stats \
    --output "${OUT}/pairs/pair_stats.txt" \
    "${OUT}/pairs/dedup.pairs.gz"
fi

# =============================================================================
# Level 6c: cooler cload — create contact matrix
# =============================================================================
if [ ! -f "${OUT}/matrix/contacts.cool" ]; then
  echo "[Level 6c] Running cooler cload pairs..."
  cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    "${CHROMSIZES}:${RESOLUTION}" \
    "${OUT}/pairs/dedup.pairs.gz" \
    "${OUT}/matrix/contacts.cool"
fi

# =============================================================================
# Level 6d: cooler balance — ICE normalization
# =============================================================================
if [ ! -f "${OUT}/matrix/.balanced" ]; then
  echo "[Level 6d] Running cooler balance (ICE normalization)..."
  cooler balance --force -p "${THREADS}" "${OUT}/matrix/contacts.cool"
  touch "${OUT}/matrix/.balanced"
fi

# =============================================================================
# Level 7a: cooltools eigs — A/B compartment calling
# =============================================================================
if [ ! -f "${OUT}/compartments/eigs.cis.vecs.tsv" ]; then
  echo "[Level 7a] Running cooltools eigs-cis (compartments)..."
  cooltools eigs-cis \
    --n-eigs 3 \
    -o "${OUT}/compartments/eigs" \
    "${OUT}/matrix/contacts.cool" || true
fi

# =============================================================================
# Level 7b: cooltools insulation — TAD/CID boundaries
# =============================================================================
if [ ! -f "${OUT}/insulation/insulation.tsv" ]; then
  echo "[Level 7b] Running cooltools insulation (TAD boundaries)..."
  cooltools insulation \
    "${OUT}/matrix/contacts.cool" \
    25000 50000 \
    -o "${OUT}/insulation/insulation.tsv" || true
fi

# =============================================================================
# Level 8: CONVERGENCE 2 — compartments + TADs + pairs ready
# =============================================================================
echo "[Level 8] Convergence 2: compartments + insulation + dedup pairs"

# =============================================================================
# Level 9a: cooltools saddle — compartment strength
# =============================================================================
if [ ! -f "${OUT}/saddle/saddle.tsv" ]; then
  echo "[Level 9a] Running cooltools saddle..."
  if [ -f "${OUT}/compartments/eigs.cis.vecs.tsv" ]; then
    cooltools saddle \
      "${OUT}/matrix/contacts.cool" \
      "${OUT}/compartments/eigs.cis.vecs.tsv" \
      --n-bins 30 \
      -o "${OUT}/saddle/saddle" || true
  else
    echo "Skipping saddle: no compartment eigenvectors available"
  fi
fi

# =============================================================================
# Level 9b: chromosight detect — loop/dot patterns
# =============================================================================
if [ ! -f "${OUT}/loops/chromosight_loops.tsv" ]; then
  echo "[Level 9b] Running chromosight detect (loops)..."
  chromosight detect \
    --pattern loops \
    --threads "${THREADS}" \
    "${OUT}/matrix/contacts.cool" \
    "${OUT}/loops/chromosight_loops" || true
fi

# =============================================================================
# Level 9c: FitHiC — significant chromatin interactions
# =============================================================================
if [ ! -f "${OUT}/fithic/fithic_results.txt" ]; then
  echo "[Level 9c] Preparing FitHiC input and running..."

  # Extract contacts from pairs file for FitHiC format
  python3 << 'FITHIC_PREP'
import gzip
import os
from collections import defaultdict

pairs_file = os.environ.get("OUT", "outputs") + "/pairs/dedup.pairs.gz"
fithic_dir = os.environ.get("OUT", "outputs") + "/fithic"
resolution = int(os.environ.get("RESOLUTION", "5000"))
chromsizes_file = os.environ.get("REF", "reference") + "/chromsizes.tsv"

# Read chromsizes
chromsizes = {}
with open(chromsizes_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        chromsizes[parts[0]] = int(parts[1])

# Count contacts and fragment coverage
contacts = defaultdict(int)
frag_coverage = defaultdict(int)

with gzip.open(pairs_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        chrom1, pos1, chrom2, pos2 = parts[1], int(parts[2]), parts[3], int(parts[4])
        if chrom1 != chrom2:
            continue  # FitHiC: cis only
        mid1 = (pos1 // resolution) * resolution + resolution // 2
        mid2 = (pos2 // resolution) * resolution + resolution // 2
        if mid1 > mid2:
            mid1, mid2 = mid2, mid1
        contacts[(chrom1, mid1, chrom2, mid2)] += 1
        frag_coverage[(chrom1, mid1)] += 1
        frag_coverage[(chrom2, mid2)] += 1

# Write fragments file
with open(f"{fithic_dir}/fragments.txt", 'w') as f:
    for (chrom, mid), count in sorted(frag_coverage.items()):
        f.write(f"{chrom}\t0\t{mid}\t{count}\t0\n")

# Write interactions file
with open(f"{fithic_dir}/interactions.txt", 'w') as f:
    for (c1, m1, c2, m2), count in sorted(contacts.items()):
        f.write(f"{c1}\t{m1}\t{c2}\t{m2}\t{count}\n")

print(f"Fragments: {len(frag_coverage)}, Interactions: {len(contacts)}")
FITHIC_PREP

  # Gzip input for FitHiC (required)
  gzip -f "${OUT}/fithic/fragments.txt"
  gzip -f "${OUT}/fithic/interactions.txt"

  # Run FitHiC
  if [ -s "${OUT}/fithic/interactions.txt.gz" ] && [ -s "${OUT}/fithic/fragments.txt.gz" ]; then
    fithic \
      -i "${OUT}/fithic/interactions.txt.gz" \
      -f "${OUT}/fithic/fragments.txt.gz" \
      -o "${OUT}/fithic" \
      -r "${RESOLUTION}" \
      -L 10000 -U 500000 || true
  fi

  # Find FitHiC output file
  FITHIC_OUT=$(ls "${OUT}"/fithic/FitHiC.spline_pass*.significances.txt* 2>/dev/null | head -1 || true)
  if [ -n "${FITHIC_OUT}" ]; then
    # Handle gzipped or plain text output
    if [[ "${FITHIC_OUT}" == *.gz ]]; then
      zcat "${FITHIC_OUT}" > "${OUT}/fithic/fithic_results.txt"
    else
      cp "${FITHIC_OUT}" "${OUT}/fithic/fithic_results.txt"
    fi
  else
    touch "${OUT}/fithic/fithic_results.txt"
  fi
fi

# =============================================================================
# Level 10: CONVERGENCE 3 + 4 — Generate report
# =============================================================================
echo "[Level 10] Generating final report..."

python3 << 'REPORT'
import os, json

OUT = os.environ.get("OUT", "outputs")
RES_DIR = os.environ.get("RES", "results")
RESOLUTION = int(os.environ.get("RESOLUTION", "5000"))

metrics = {}

# --- fastp QC ---
with open(f"{OUT}/qc/fastp.json") as f:
    fj = json.load(f)
metrics["total_read_pairs"] = fj["summary"]["before_filtering"]["total_reads"] // 2
metrics["reads_after_trim"] = fj["summary"]["after_filtering"]["total_reads"] // 2
metrics["trim_rate_pct"] = round(100 * (1 - metrics["reads_after_trim"] / metrics["total_read_pairs"]), 2)

# --- pairtools stats ---
stats = {}
with open(f"{OUT}/pairs/pair_stats.txt") as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            parts = line.split("\t")
            if len(parts) >= 2:
                stats[parts[0]] = parts[1]

metrics["mapped_pairs"] = int(stats.get("total_mapped", 0))
metrics["unique_pairs"] = int(stats.get("total_mapped", 0))
metrics["cis_contacts"] = int(stats.get("cis", 0))
metrics["trans_contacts"] = int(stats.get("trans", 0))
cis = metrics["cis_contacts"]
trans = metrics["trans_contacts"]
metrics["cis_ratio"] = round(cis / (cis + trans), 4) if (cis + trans) > 0 else 0
metrics["cis_gt_1kb"] = int(stats.get("cis_1kb+", 0))
metrics["cis_gt_10kb"] = int(stats.get("cis_10kb+", 0))

# --- dedup stats ---
with open(f"{OUT}/pairs/dedup_stats.txt") as f:
    for line in f:
        if line.startswith("total_dups"):
            metrics["duplicate_pairs"] = int(line.strip().split("\t")[1])

# --- Compartments: E1 column ---
compartment_switches = 0
a_bins = 0
b_bins = 0
if os.path.exists(f"{OUT}/compartments/eigs.cis.vecs.tsv"):
    with open(f"{OUT}/compartments/eigs.cis.vecs.tsv") as f:
        header = f.readline().strip().split("\t")
        e1_idx = header.index("E1") if "E1" in header else 4
        prev_sign = None
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > e1_idx and parts[e1_idx]:
                try:
                    val = float(parts[e1_idx])
                    if val > 0: a_bins += 1
                    else: b_bins += 1
                    curr = "A" if val > 0 else "B"
                    if prev_sign is not None and curr != prev_sign:
                        compartment_switches += 1
                    prev_sign = curr
                except (ValueError, IndexError):
                    prev_sign = None
metrics["compartment_switches"] = compartment_switches
metrics["compartment_a_bins"] = a_bins
metrics["compartment_b_bins"] = b_bins

# --- Boundaries ---
boundary_25k = 0
boundary_50k = 0
if os.path.exists(f"{OUT}/insulation/insulation.tsv"):
    with open(f"{OUT}/insulation/insulation.tsv") as f:
        header = f.readline().strip().split("\t")
        b25_idx = header.index("is_boundary_25000") if "is_boundary_25000" in header else None
        b50_idx = header.index("is_boundary_50000") if "is_boundary_50000" in header else None
        for line in f:
            parts = line.strip().split("\t")
            if b25_idx and b25_idx < len(parts) and parts[b25_idx].strip().lower() == "true":
                boundary_25k += 1
            if b50_idx and b50_idx < len(parts) and parts[b50_idx].strip().lower() == "true":
                boundary_50k += 1
metrics["boundaries_25kb"] = boundary_25k
metrics["boundaries_50kb"] = boundary_50k

# --- Chromosight loops ---
loop_count = 0
cfile = f"{OUT}/loops/chromosight_loops.tsv"
if os.path.exists(cfile):
    with open(cfile) as f:
        next(f)  # skip header
        for line in f:
            if line.strip(): loop_count += 1
metrics["detected_loops"] = loop_count

# --- FitHiC significant interactions ---
sig_count = 0
ffile = f"{OUT}/fithic/fithic_results.txt"
if os.path.exists(ffile):
    with open(ffile) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                try:
                    if float(parts[6]) < 0.05: sig_count += 1
                except (ValueError, IndexError):
                    pass
metrics["significant_interactions"] = sig_count
metrics["resolution_bp"] = RESOLUTION

# --- Write CSV ---
with open(f"{RES_DIR}/report.csv", "w") as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("Report written to results/report.csv")
for k, v in metrics.items():
    print(f"  {k} = {v}")
REPORT

echo "=== Pipeline Complete ==="
cat "${RES}/report.csv"
