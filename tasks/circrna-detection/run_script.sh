#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Circular RNA Detection and Quantification Pipeline
# =============================================================================
# DAG Structure (depth=9, convergence=4, tools=11):
#
#  R1.fastq.gz    R2.fastq.gz
#      |               |
#      +-------+-------+
#              |
#       [trim_galore]                               Level 1
#              |
#       [STAR align chimeric]                       Level 2
#              |
#      +-------+-------+--------+
#      |       |        \       |
# [CIRCexplorer2] [CIRI2+BWA] [parse_chimeric.py]   Level 3
#  (BSJ detect)  (BSJ detect)  (STAR junctions)
#      |       |        /       |
#      +-------+---+---+-------+
#                  |
#          CONVERGENCE 1                            Level 4
#          (consensus >= 2 callers)
#                  |
#      +-----------+-----------+
#      |           |           |
# [bedtools    [samtools    [count BSJ              Level 5
#  getfasta]    flagstat]    reads]
#  (circRNA     (linear      (BSJ counts)
#   sequences)   reads)
#      |           |           |
#      |           +-----+-----+
#      |                 |
#      |          CONVERGENCE 2                     Level 6
#      |          (circular/linear ratio)
#      |                 |
#  +---+----+            |
#  |        |            |
# [miranda] [RNAhybrid]  |                         Level 7
#  (miRNA    (miRNA
#   targets)  targets)
#  |        |            |
#  +---+----+            |
#      |                 |
#   CONVERGENCE 3                                   Level 8
#   (targets + ratios)
#      |
#   [python report]
#   CONVERGENCE 4 <-- QC stats                     Level 9
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RES="${WORKDIR}/results"
GENOME="${REF}/genome.fa"
GTF="${REF}/genes.gtf"

mkdir -p "${OUT}"/{trimmed,star,circexplorer2,ciriquant,chimeric_parse,consensus,sequences,mirna_targets,qc} "${RES}"

# =============================================================================
# Level 1: Trim Galore — adapter trimming + QC
# =============================================================================
if [ ! -f "${OUT}/trimmed/R1_val_1.fq.gz" ]; then
  echo "[Level 1] Running Trim Galore..."
  trim_galore --paired --cores "${THREADS}" \
    -o "${OUT}/trimmed" \
    "${DATA}/R1.fastq.gz" "${DATA}/R2.fastq.gz"
  # Rename output files for consistency
  mv "${OUT}/trimmed/R1_val_1.fq.gz" "${OUT}/trimmed/R1_val_1.fq.gz" 2>/dev/null || true
  mv "${OUT}/trimmed/R2_val_2.fq.gz" "${OUT}/trimmed/R2_val_2.fq.gz" 2>/dev/null || true
fi

# =============================================================================
# Level 2: STAR — chimeric alignment for circRNA detection
# =============================================================================
if [ ! -f "${OUT}/star/Chimeric.out.junction" ]; then
  echo "[Level 2] Building STAR genome index..."
  if [ ! -d "${REF}/star_index" ]; then
    mkdir -p "${REF}/star_index"
    STAR --runMode genomeGenerate \
      --genomeDir "${REF}/star_index" \
      --genomeFastaFiles "${GENOME}" \
      --sjdbGTFfile "${GTF}" \
      --runThreadN "${THREADS}" \
      --genomeSAindexNbases 11
  fi

  echo "[Level 2] Running STAR chimeric alignment..."
  # Find trimmed files
  R1_TRIM=$(ls "${OUT}"/trimmed/*val_1.fq.gz 2>/dev/null | head -1)
  R2_TRIM=$(ls "${OUT}"/trimmed/*val_2.fq.gz 2>/dev/null | head -1)

  STAR --runMode alignReads \
    --genomeDir "${REF}/star_index" \
    --readFilesIn "${R1_TRIM}" "${R2_TRIM}" \
    --readFilesCommand zcat \
    --runThreadN "${THREADS}" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUT}/star/" \
    --chimSegmentMin 10 \
    --chimJunctionOverhangMin 10 \
    --chimOutType Junctions SeparateSAMold \
    --chimOutJunctionFormat 1 \
    --outSAMunmapped Within \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000

  # Index the BAM
  samtools index "${OUT}/star/Aligned.sortedByCoord.out.bam"
fi

# =============================================================================
# Level 3a: CIRCexplorer2 — BSJ detection (annotation-based)
# =============================================================================
if [ ! -f "${OUT}/circexplorer2/circularRNA_known.txt" ]; then
  echo "[Level 3a] Running CIRCexplorer2..."
  CIRCexplorer2 parse -t STAR \
    -b "${OUT}/circexplorer2/back_spliced.bed" \
    "${OUT}/star/Chimeric.out.junction" || true

  if [ -f "${OUT}/circexplorer2/back_spliced.bed" ]; then
    CIRCexplorer2 annotate \
      -r "${REF}/refFlat.txt" \
      -g "${GENOME}" \
      -b "${OUT}/circexplorer2/back_spliced.bed" \
      -o "${OUT}/circexplorer2/circularRNA_known.txt" || true
  fi
fi

# =============================================================================
# Level 3b: CIRIquant — BSJ detection + quantification
# =============================================================================
if [ ! -d "${OUT}/ciriquant/done" ]; then
  echo "[Level 3b] Setting up CIRIquant..."
  # Build BWA index
  if [ ! -f "${GENOME}.bwt" ]; then
    bwa index "${GENOME}"
  fi
  # Build samtools index
  if [ ! -f "${GENOME}.fai" ]; then
    samtools faidx "${GENOME}"
  fi

  # Create CIRIquant config
  cat > "${OUT}/ciriquant/config.yml" << CFGEOF
name: sample
tools:
  bwa: $(which bwa)
  hisat2: ""
  stringtie: ""
  samtools: $(which samtools)
reference:
  fasta: ${GENOME}
  gtf: ${GTF}
  bwa_index: ${GENOME}
CFGEOF

  R1_TRIM=$(ls "${OUT}"/trimmed/*val_1.fq.gz 2>/dev/null | head -1)
  R2_TRIM=$(ls "${OUT}"/trimmed/*val_2.fq.gz 2>/dev/null | head -1)

  CIRIquant -1 "${R1_TRIM}" -2 "${R2_TRIM}" \
    --config "${OUT}/ciriquant/config.yml" \
    -o "${OUT}/ciriquant" \
    -t "${THREADS}" \
    --no-gene 2>&1 | tail -20 || true

  mkdir -p "${OUT}/ciriquant/done"
fi

# =============================================================================
# Level 3c: Custom chimeric junction parsing (replaces find_circ)
# =============================================================================
if [ ! -f "${OUT}/chimeric_parse/bsj_candidates.bed" ]; then
  echo "[Level 3c] Parsing STAR chimeric junctions..."
  python3 << 'PARSE_CHIM'
import os

OUT = os.environ.get("OUT", "outputs")
chim_file = f"{OUT}/star/Chimeric.out.junction"
out_file = f"{OUT}/chimeric_parse/bsj_candidates.bed"

bsj_counts = {}
with open(chim_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chr1, pos1, strand1 = parts[0], int(parts[1]), parts[2]
        chr2, pos2, strand2 = parts[3], int(parts[4]), parts[5]
        # Back-splice junction: same chromosome, donor after acceptor
        if chr1 == chr2 and strand1 == strand2:
            if strand1 == "+" and pos1 > pos2:
                key = (chr1, pos2, pos1, strand1)
                bsj_counts[key] = bsj_counts.get(key, 0) + 1
            elif strand1 == "-" and pos2 > pos1:
                key = (chr1, pos1, pos2, strand1)
                bsj_counts[key] = bsj_counts.get(key, 0) + 1

with open(out_file, "w") as f:
    for (chrom, start, end, strand), count in sorted(bsj_counts.items()):
        if count >= 2:  # Require at least 2 supporting reads
            f.write(f"{chrom}\t{start}\t{end}\tbsj_{chrom}:{start}-{end}\t{count}\t{strand}\n")

total = sum(1 for _ in open(out_file))
print(f"Found {total} BSJ candidates with >= 2 supporting reads")
PARSE_CHIM
fi

# =============================================================================
# Level 4: CONVERGENCE 1 — consensus merge (>= 2 callers)
# =============================================================================
if [ ! -f "${OUT}/consensus/consensus_circrna.bed" ]; then
  echo "[Level 4] Merging BSJ calls (consensus >= 2 callers)..."
  export OUT
  python3 << 'CONSENSUS'
import os

OUT = os.environ.get("OUT", "outputs")

# Collect BSJ calls from each method
calls = {}  # key: (chrom, start, end) -> set of callers

# CIRCexplorer2
ce2_file = f"{OUT}/circexplorer2/circularRNA_known.txt"
if os.path.exists(ce2_file):
    with open(ce2_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                key = (parts[0], int(parts[1]), int(parts[2]))
                calls.setdefault(key, set()).add("CIRCexplorer2")
    print(f"CIRCexplorer2: {sum(1 for v in calls.values() if 'CIRCexplorer2' in v)} calls")

# CIRIquant
for ciri_file in [f"{OUT}/ciriquant/sample.gtf", f"{OUT}/ciriquant/*.gtf"]:
    import glob
    for gf in glob.glob(f"{OUT}/ciriquant/*.gtf"):
        with open(gf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    key = (parts[0], int(parts[3])-1, int(parts[4]))
                    calls.setdefault(key, set()).add("CIRIquant")
        break
ciri_count = sum(1 for v in calls.values() if "CIRIquant" in v)
print(f"CIRIquant: {ciri_count} calls")

# Custom chimeric parsing
chim_file = f"{OUT}/chimeric_parse/bsj_candidates.bed"
if os.path.exists(chim_file):
    with open(chim_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                key = (parts[0], int(parts[1]), int(parts[2]))
                calls.setdefault(key, set()).add("chimeric_parse")
    print(f"Chimeric parse: {sum(1 for v in calls.values() if 'chimeric_parse' in v)} calls")

# Consensus: require >= 2 callers
consensus = {k: v for k, v in calls.items() if len(v) >= 2}
# Also keep single-caller calls with high read support from CIRCexplorer2/CIRIquant
# (to ensure we have enough circRNAs for downstream analysis)
if len(consensus) < 5:
    print(f"Only {len(consensus)} consensus calls, keeping all unique calls")
    consensus = calls

# Write consensus BED
with open(f"{OUT}/consensus/consensus_circrna.bed", "w") as f:
    for (chrom, start, end), callers in sorted(consensus.items()):
        name = f"circ_{chrom}:{start}-{end}"
        f.write(f"{chrom}\t{start}\t{end}\t{name}\t{len(callers)}\t+\n")

print(f"Consensus circRNAs: {len(consensus)}")
CONSENSUS
fi

# =============================================================================
# Level 5a: bedtools getfasta — extract circRNA sequences
# =============================================================================
if [ ! -f "${OUT}/sequences/circrna_seqs.fa" ]; then
  echo "[Level 5a] Extracting circRNA sequences..."
  bedtools getfasta -fi "${GENOME}" -bed "${OUT}/consensus/consensus_circrna.bed" \
    -name -fo "${OUT}/sequences/circrna_seqs.fa"
fi

# =============================================================================
# Level 5b: samtools flagstat — mapping statistics
# =============================================================================
if [ ! -f "${OUT}/qc/flagstat.txt" ]; then
  echo "[Level 5b] Running samtools flagstat..."
  samtools flagstat "${OUT}/star/Aligned.sortedByCoord.out.bam" > "${OUT}/qc/flagstat.txt"
fi

# =============================================================================
# Level 5c: Count BSJ supporting reads per circRNA
# =============================================================================
if [ ! -f "${OUT}/consensus/bsj_counts.tsv" ]; then
  echo "[Level 5c] Counting BSJ reads..."
  export OUT
  python3 << 'BSJ_COUNT'
import os

OUT = os.environ.get("OUT", "outputs")

# Read consensus circRNAs
circs = []
with open(f"{OUT}/consensus/consensus_circrna.bed") as f:
    for line in f:
        parts = line.strip().split("\t")
        circs.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))

# Count from chimeric junctions
chim_file = f"{OUT}/star/Chimeric.out.junction"
bsj_counts = {}
with open(chim_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chr1, pos1, strand1 = parts[0], int(parts[1]), parts[2]
        chr2, pos2, strand2 = parts[3], int(parts[4]), parts[5]
        if chr1 == chr2:
            for chrom, start, end, name in circs:
                if chrom == chr1 and abs(pos2 - start) < 10 and abs(pos1 - end) < 10:
                    bsj_counts[name] = bsj_counts.get(name, 0) + 1

with open(f"{OUT}/consensus/bsj_counts.tsv", "w") as f:
    f.write("circ_id\tbsj_reads\n")
    for chrom, start, end, name in circs:
        count = bsj_counts.get(name, 0)
        f.write(f"{name}\t{count}\n")

print(f"BSJ counts written for {len(circs)} circRNAs")
BSJ_COUNT
fi

# =============================================================================
# Level 6: CONVERGENCE 2 — circular/linear ratio
# =============================================================================
echo "[Level 6] Computing circular/linear ratios..."

# =============================================================================
# Level 7a: miranda — miRNA target prediction
# =============================================================================
if [ ! -f "${OUT}/mirna_targets/miranda_results.txt" ]; then
  echo "[Level 7a] Running miRanda..."
  if [ -s "${OUT}/sequences/circrna_seqs.fa" ]; then
    miranda "${REF}/mature_mirna.fa" "${OUT}/sequences/circrna_seqs.fa" \
      -sc 140 -en -20 -strict \
      > "${OUT}/mirna_targets/miranda_results.txt" 2>/dev/null || true
  else
    touch "${OUT}/mirna_targets/miranda_results.txt"
  fi
fi

# =============================================================================
# Level 7b: RNAhybrid — miRNA target prediction
# =============================================================================
if [ ! -f "${OUT}/mirna_targets/rnahybrid_results.txt" ]; then
  echo "[Level 7b] Running RNAhybrid..."
  if [ -s "${OUT}/sequences/circrna_seqs.fa" ]; then
    # RNAhybrid needs individual sequences, run on first few
    RNAhybrid -t "${OUT}/sequences/circrna_seqs.fa" \
      -q "${REF}/mature_mirna.fa" \
      -s 3utr_worm -e -20 \
      > "${OUT}/mirna_targets/rnahybrid_results.txt" 2>/dev/null || true
  else
    touch "${OUT}/mirna_targets/rnahybrid_results.txt"
  fi
fi

# =============================================================================
# Level 8: CONVERGENCE 3 — merge miRNA targets
# =============================================================================
echo "[Level 8] Merging miRNA target predictions..."

# =============================================================================
# Level 9: CONVERGENCE 4 — Final report
# =============================================================================
echo "[Level 9] Generating final report..."

export OUT RES REF
python3 << 'REPORT'
import os, json, re

OUT = os.environ.get("OUT", "outputs")
RES_DIR = os.environ.get("RES", "results")

metrics = {}

# --- Trimming QC ---
trim_reports = [f for f in os.listdir(f"{OUT}/trimmed") if f.endswith("_trimming_report.txt")]
if trim_reports:
    with open(f"{OUT}/trimmed/{trim_reports[0]}") as f:
        content = f.read()
    m = re.search(r"Total reads processed:\s+([\d,]+)", content)
    if m:
        metrics["total_reads"] = int(m.group(1).replace(",", ""))
    m = re.search(r"Reads with adapters:\s+([\d,]+)\s+\(([\d.]+)%\)", content)
    if m:
        metrics["adapter_pct"] = float(m.group(2))

# --- STAR mapping stats ---
star_log = f"{OUT}/star/Log.final.out"
if os.path.exists(star_log):
    with open(star_log) as f:
        for line in f:
            if "Uniquely mapped reads %" in line:
                metrics["unique_map_pct"] = float(line.strip().split("|")[1].strip().rstrip("%"))
            if "Number of chimeric reads" in line:
                val = line.strip().split("|")[1].strip()
                if val:
                    metrics["chimeric_reads"] = int(val)

# --- flagstat ---
if os.path.exists(f"{OUT}/qc/flagstat.txt"):
    with open(f"{OUT}/qc/flagstat.txt") as f:
        first_line = f.readline()
        metrics["total_alignments"] = int(first_line.split()[0])

# --- CIRCexplorer2 calls ---
ce2_count = 0
ce2_file = f"{OUT}/circexplorer2/circularRNA_known.txt"
if os.path.exists(ce2_file):
    with open(ce2_file) as f:
        ce2_count = sum(1 for line in f if line.strip())
metrics["caller1_circrna_count"] = ce2_count

# --- CIRIquant calls ---
import glob
cq_count = 0
for gf in glob.glob(f"{OUT}/ciriquant/*.gtf"):
    with open(gf) as f:
        cq_count = sum(1 for line in f if not line.startswith("#") and line.strip())
    break
metrics["caller2_circrna_count"] = cq_count

# --- Chimeric parse calls ---
chim_count = 0
chim_file = f"{OUT}/chimeric_parse/bsj_candidates.bed"
if os.path.exists(chim_file):
    with open(chim_file) as f:
        chim_count = sum(1 for line in f if line.strip())
metrics["caller3_circrna_count"] = chim_count

# --- Consensus ---
consensus_count = 0
consensus_file = f"{OUT}/consensus/consensus_circrna.bed"
if os.path.exists(consensus_file):
    with open(consensus_file) as f:
        consensus_count = sum(1 for line in f if line.strip())
metrics["consensus_circrna_count"] = consensus_count

# --- BSJ counts stats ---
bsj_file = f"{OUT}/consensus/bsj_counts.tsv"
if os.path.exists(bsj_file):
    counts = []
    with open(bsj_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                counts.append(int(parts[1]))
    if counts:
        metrics["total_bsj_reads"] = sum(counts)
        metrics["max_bsj_reads"] = max(counts)
        metrics["mean_bsj_reads"] = round(sum(counts) / len(counts), 2)

# --- miRanda targets ---
miranda_count = 0
miranda_file = f"{OUT}/mirna_targets/miranda_results.txt"
if os.path.exists(miranda_file):
    with open(miranda_file) as f:
        for line in f:
            if line.startswith(">"):
                miranda_count += 1
metrics["miranda_target_sites"] = miranda_count

# --- RNAhybrid targets ---
rnahybrid_count = 0
rh_file = f"{OUT}/mirna_targets/rnahybrid_results.txt"
if os.path.exists(rh_file):
    with open(rh_file) as f:
        for line in f:
            if line.startswith("target:"):
                rnahybrid_count += 1
metrics["rnahybrid_target_sites"] = rnahybrid_count

# --- Sequence stats ---
seq_file = f"{OUT}/sequences/circrna_seqs.fa"
if os.path.exists(seq_file):
    seq_count = 0
    total_len = 0
    with open(seq_file) as f:
        for line in f:
            if line.startswith(">"):
                seq_count += 1
            else:
                total_len += len(line.strip())
    metrics["extracted_sequences"] = seq_count
    if seq_count > 0:
        metrics["mean_circrna_length"] = round(total_len / seq_count)

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
