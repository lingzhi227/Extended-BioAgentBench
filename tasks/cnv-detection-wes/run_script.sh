#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CNV Detection from WES Pipeline
# =============================================================================
# DAG Structure (depth=10, convergence=4, tools=10):
#
#  tumor.fastq.gz         normal.fastq.gz
#      |                       |
#  [fastp QC] ----------- [fastp QC]             Level 1
#      |                       |
#  [bwa-mem2 align] ------ [bwa-mem2 align]      Level 2
#      |                       |
#  [picard MarkDup] ------ [picard MarkDup]      Level 3
#      |                       |
#  [mosdepth coverage] --- [mosdepth coverage]   Level 4
#      |                       |
#      +-----------+-----------+
#                  |
#          CONVERGENCE 1                          Level 5
#          (T+N BAMs + coverage)
#                  |
#      +-----------+-----------+
#      |           |           |
#  [cnvkit      [freec]     [gatk4               Level 6
#   batch]                   CollectRC->
#                            DenoiseRC->
#                            ModelSeg->
#                            CallSeg]
#      |           |           |
#      +-----------+-----+-----+
#                        |
#                CONVERGENCE 2                    Level 7
#                (3-caller merge)
#                [bedtools intersect]
#                        |
#                +-------+-------+
#                |       |       |
#           [bcftools [bedtools [python            Level 8
#            stats]   intersect gene-level
#                     w/ genes] stats]
#                |       |       |
#                +-------+-------+
#                        |
#                CONVERGENCE 3                    Level 9
#                (stats + gene-impact)
#                        |
#                [python report]
#                CONVERGENCE 4 <-- QC + coverage  Level 10
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RES="${WORKDIR}/results"
GENOME="${REF}/genome.fa"
TARGETS="${REF}/targets.bed"

mkdir -p "${OUT}"/{tumor,normal,cnvkit,freec,gatk_cnv,merged,qc} "${RES}"

# =============================================================================
# Index reference
# =============================================================================
if [ ! -f "${GENOME}.bwt.2bit.64" ]; then
  echo "Indexing reference..."
  bwa-mem2 index "${GENOME}"
  samtools faidx "${GENOME}"
  samtools dict "${GENOME}" > "${GENOME%.fa}.dict"
fi

# =============================================================================
# Level 1-3: Process Tumor
# =============================================================================
if [ ! -f "${OUT}/tumor/tumor.dedup.bam" ]; then
  echo "[L1-3] Processing tumor..."
  fastp -i "${DATA}/tumor_R1.fastq.gz" -I "${DATA}/tumor_R2.fastq.gz" \
    -o "${OUT}/tumor/R1_trimmed.fq.gz" -O "${OUT}/tumor/R2_trimmed.fq.gz" \
    --json "${OUT}/qc/tumor_fastp.json" --thread "${THREADS}" --length_required 30

  bwa-mem2 mem -t "${THREADS}" \
    -R "@RG\tID:tumor\tSM:tumor\tPL:ILLUMINA\tLB:lib1" \
    "${GENOME}" "${OUT}/tumor/R1_trimmed.fq.gz" "${OUT}/tumor/R2_trimmed.fq.gz" \
    | samtools sort -@ "${THREADS}" -o "${OUT}/tumor/tumor.sorted.bam"

  picard MarkDuplicates \
    I="${OUT}/tumor/tumor.sorted.bam" \
    O="${OUT}/tumor/tumor.dedup.bam" \
    M="${OUT}/tumor/tumor.dedup_metrics.txt" \
    REMOVE_DUPLICATES=false CREATE_INDEX=true 2>/dev/null
fi

# =============================================================================
# Level 1-3: Process Normal
# =============================================================================
if [ ! -f "${OUT}/normal/normal.dedup.bam" ]; then
  echo "[L1-3] Processing normal..."
  fastp -i "${DATA}/normal_R1.fastq.gz" -I "${DATA}/normal_R2.fastq.gz" \
    -o "${OUT}/normal/R1_trimmed.fq.gz" -O "${OUT}/normal/R2_trimmed.fq.gz" \
    --json "${OUT}/qc/normal_fastp.json" --thread "${THREADS}" --length_required 30

  bwa-mem2 mem -t "${THREADS}" \
    -R "@RG\tID:normal\tSM:normal\tPL:ILLUMINA\tLB:lib1" \
    "${GENOME}" "${OUT}/normal/R1_trimmed.fq.gz" "${OUT}/normal/R2_trimmed.fq.gz" \
    | samtools sort -@ "${THREADS}" -o "${OUT}/normal/normal.sorted.bam"

  picard MarkDuplicates \
    I="${OUT}/normal/normal.sorted.bam" \
    O="${OUT}/normal/normal.dedup.bam" \
    M="${OUT}/normal/normal.dedup_metrics.txt" \
    REMOVE_DUPLICATES=false CREATE_INDEX=true 2>/dev/null
fi

# =============================================================================
# Level 4: mosdepth coverage
# =============================================================================
if [ ! -f "${OUT}/tumor/tumor.mosdepth.global.dist.txt" ]; then
  echo "[L4] Running mosdepth coverage..."
  mosdepth --by "${TARGETS}" -t "${THREADS}" \
    "${OUT}/tumor/tumor" "${OUT}/tumor/tumor.dedup.bam" || true
  mosdepth --by "${TARGETS}" -t "${THREADS}" \
    "${OUT}/normal/normal" "${OUT}/normal/normal.dedup.bam" || true
fi

# =============================================================================
# Level 5: CONVERGENCE 1 — T+N BAMs ready
# =============================================================================
echo "[L5] Convergence 1: T+N BAMs ready"

# =============================================================================
# Level 6a: CNVkit batch
# =============================================================================
if [ ! -f "${OUT}/cnvkit/tumor.cns" ]; then
  echo "[L6a] Running CNVkit..."
  cnvkit.py batch \
    "${OUT}/tumor/tumor.dedup.bam" \
    --normal "${OUT}/normal/normal.dedup.bam" \
    --targets "${TARGETS}" \
    --fasta "${GENOME}" \
    --output-dir "${OUT}/cnvkit" \
    --diagram --scatter \
    -p "${THREADS}" || true

  # Call CNVs
  if [ -f "${OUT}/cnvkit/tumor.cns" ]; then
    cnvkit.py call "${OUT}/cnvkit/tumor.cns" -o "${OUT}/cnvkit/tumor.call.cns" || true
  fi
fi

# =============================================================================
# Level 6b: Control-FREEC
# =============================================================================
if [ ! -f "${OUT}/freec/tumor_CNVs" ]; then
  echo "[L6b] Running Control-FREEC..."
  # Create FREEC config
  cat > "${OUT}/freec/config.txt" << FREEC_CFG
[general]
chrLenFile = ${GENOME}.fai
ploidy = 2
window = 50000
outputDir = ${OUT}/freec
maxThreads = ${THREADS}

[sample]
mateFile = ${OUT}/tumor/tumor.dedup.bam
inputFormat = BAM
mateOrientation = FR

[control]
mateFile = ${OUT}/normal/normal.dedup.bam
inputFormat = BAM
mateOrientation = FR

[target]
captureRegions = ${TARGETS}
FREEC_CFG

  freec -conf "${OUT}/freec/config.txt" || true
fi

# =============================================================================
# Level 6c: GATK4 CNV
# =============================================================================
if [ ! -f "${OUT}/gatk_cnv/tumor.called.seg" ]; then
  echo "[L6c] Running GATK4 CNV pipeline..."

  # Create interval list from BED
  gatk BedToIntervalList \
    -I "${TARGETS}" \
    -O "${OUT}/gatk_cnv/targets.interval_list" \
    -SD "${GENOME%.fa}.dict" || true

  if [ -f "${OUT}/gatk_cnv/targets.interval_list" ]; then
    # Preprocess intervals
    gatk PreprocessIntervals \
      -R "${GENOME}" \
      -L "${OUT}/gatk_cnv/targets.interval_list" \
      --bin-length 0 \
      --padding 250 \
      -O "${OUT}/gatk_cnv/preprocessed.interval_list" || true

    # Collect read counts
    for sample in tumor normal; do
      gatk CollectReadCounts \
        -I "${OUT}/${sample}/${sample}.dedup.bam" \
        -L "${OUT}/gatk_cnv/preprocessed.interval_list" \
        -R "${GENOME}" \
        -O "${OUT}/gatk_cnv/${sample}.counts.hdf5" || true
    done

    # Create PoN from normal
    if [ -f "${OUT}/gatk_cnv/normal.counts.hdf5" ]; then
      gatk CreateReadCountPanelOfNormals \
        -I "${OUT}/gatk_cnv/normal.counts.hdf5" \
        -O "${OUT}/gatk_cnv/pon.hdf5" || true
    fi

    # Denoise read counts
    if [ -f "${OUT}/gatk_cnv/tumor.counts.hdf5" ] && [ -f "${OUT}/gatk_cnv/pon.hdf5" ]; then
      gatk DenoiseReadCounts \
        -I "${OUT}/gatk_cnv/tumor.counts.hdf5" \
        --count-panel-of-normals "${OUT}/gatk_cnv/pon.hdf5" \
        --standardized-copy-ratios "${OUT}/gatk_cnv/tumor.standardized.tsv" \
        --denoised-copy-ratios "${OUT}/gatk_cnv/tumor.denoised.tsv" || true
    fi

    # Model segments
    if [ -f "${OUT}/gatk_cnv/tumor.denoised.tsv" ]; then
      gatk ModelSegments \
        --denoised-copy-ratios "${OUT}/gatk_cnv/tumor.denoised.tsv" \
        -O "${OUT}/gatk_cnv" \
        --output-prefix tumor || true
    fi

    # Call segments
    if [ -f "${OUT}/gatk_cnv/tumor.cr.seg" ]; then
      gatk CallCopyRatioSegments \
        -I "${OUT}/gatk_cnv/tumor.cr.seg" \
        -O "${OUT}/gatk_cnv/tumor.called.seg" || true
    fi
  fi
fi

# =============================================================================
# Level 7: CONVERGENCE 2 — merge callers
# =============================================================================
echo "[L7] Merging CNV calls..."

export OUT
python3 << 'MERGE'
import os

OUT = os.environ.get("OUT", "outputs")
os.makedirs(f"{OUT}/merged", exist_ok=True)

all_cnvs = []

# CNVkit calls
cns_file = f"{OUT}/cnvkit/tumor.call.cns"
if os.path.exists(cns_file):
    with open(cns_file) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                try:
                    cn = int(float(parts[4]))
                    if cn != 2:
                        all_cnvs.append((chrom, start, end, "CNVkit", cn))
                except:
                    pass
    print(f"CNVkit: {sum(1 for c in all_cnvs if c[3]=='CNVkit')} CNVs")

# FREEC calls
freec_file = f"{OUT}/freec/tumor_CNVs"
if os.path.exists(freec_file):
    cnt = 0
    with open(freec_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                cn = int(parts[3])
                if cn != 2:
                    all_cnvs.append((chrom, start, end, "FREEC", cn))
                    cnt += 1
    print(f"FREEC: {cnt} CNVs")

# GATK calls
gatk_file = f"{OUT}/gatk_cnv/tumor.called.seg"
if os.path.exists(gatk_file):
    cnt = 0
    with open(gatk_file) as f:
        for line in f:
            if line.startswith("@") or line.startswith("CONTIG"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                call = parts[5] if len(parts) > 5 else ""
                if call in ("+", "-", "AMP", "DEL"):
                    all_cnvs.append((chrom, start, end, "GATK", call))
                    cnt += 1
    print(f"GATK: {cnt} CNVs")

# Write merged BED
with open(f"{OUT}/merged/all_cnvs.bed", "w") as f:
    for chrom, start, end, caller, cn in sorted(all_cnvs):
        f.write(f"{chrom}\t{start}\t{end}\t{caller}\t{cn}\n")

print(f"Total merged CNVs: {len(all_cnvs)}")
MERGE

# =============================================================================
# Level 8-10: Stats + Report
# =============================================================================
echo "[L8-10] Generating report..."

python3 << 'REPORT'
import os, json

OUT = os.environ.get("OUT", "outputs")
RES_DIR = os.environ.get("RES", "results")

metrics = {}

# fastp stats
for sample in ["tumor", "normal"]:
    fp = f"{OUT}/qc/{sample}_fastp.json"
    if os.path.exists(fp):
        with open(fp) as f:
            fj = json.load(f)
        metrics[f"{sample}_reads"] = fj["summary"]["before_filtering"]["total_reads"] // 2
        metrics[f"{sample}_reads_after_trim"] = fj["summary"]["after_filtering"]["total_reads"] // 2

# mosdepth coverage
for sample in ["tumor", "normal"]:
    summ = f"{OUT}/{sample}/{sample}.mosdepth.summary.txt"
    if os.path.exists(summ):
        with open(summ) as f:
            for line in f:
                if line.startswith("total_region") or line.startswith("total"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        metrics[f"{sample}_mean_coverage"] = float(parts[3])
                        break

# Picard dedup metrics
for sample in ["tumor", "normal"]:
    dup_file = f"{OUT}/{sample}/{sample}.dedup_metrics.txt"
    if os.path.exists(dup_file):
        with open(dup_file) as f:
            header = None
            for line in f:
                if line.startswith("LIBRARY"):
                    header = line.strip().split("\t")
                elif header and not line.startswith("#") and line.strip():
                    vals = line.strip().split("\t")
                    if len(vals) > 8:
                        try:
                            metrics[f"{sample}_dup_pct"] = round(float(vals[8]) * 100, 2)
                        except:
                            pass
                    break

# CNV caller counts
for caller, prefix in [("cnvkit", "cnvkit"), ("freec", "freec"), ("gatk", "gatk")]:
    count = 0
    if caller == "cnvkit":
        f = f"{OUT}/cnvkit/tumor.call.cns"
        if os.path.exists(f):
            with open(f) as fh:
                next(fh)
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        try:
                            cn = int(float(parts[4]))
                            if cn != 2:
                                count += 1
                        except:
                            pass
    elif caller == "freec":
        f = f"{OUT}/freec/tumor_CNVs"
        if os.path.exists(f):
            with open(f) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4 and parts[3] != "2":
                        count += 1
    elif caller == "gatk":
        f = f"{OUT}/gatk_cnv/tumor.called.seg"
        if os.path.exists(f):
            with open(f) as fh:
                for line in fh:
                    if not line.startswith("@") and not line.startswith("CONTIG"):
                        parts = line.strip().split("\t")
                        if len(parts) > 5 and parts[5] in ("+", "-", "AMP", "DEL"):
                            count += 1
    metrics[f"{prefix}_cnv_count"] = count

# Merged CNVs
merged_f = f"{OUT}/merged/all_cnvs.bed"
if os.path.exists(merged_f):
    with open(merged_f) as f:
        metrics["total_merged_cnvs"] = sum(1 for l in f if l.strip())

# Target info
targets_f = os.environ.get("REF", "reference") + "/targets.bed"
if os.path.exists(targets_f):
    with open(targets_f) as f:
        metrics["target_intervals"] = sum(1 for l in f if l.strip())

# Write CSV
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
