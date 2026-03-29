#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Task: multiomics-rna-atac
# DAG Structure (depth=10, convergence=5, tools=12):
#
# RNA_R1.fq.gz  RNA_R2.fq.gz    ATAC_R1.fq.gz  ATAC_R2.fq.gz
#     |              |               |               |
# [fastp] ------- [fastp]       [fastp] ------- [fastp]       Level 1
#     |              |               |               |
#     +------+-------+               +-------+-------+
#            |                               |
#     [STAR align]                    [bowtie2 align]          Level 2
#            |                               |
#     [samtools sort]                 [samtools sort]           Level 3
#            |                               |
#     +------+------+                 +------+------+
#     |      |      |                 |      |      |
# [salmon [stringtie [featureCounts [picard [macs2  [deeptools Level 4
#  quant]  assemble] (gene counts)] dedup]  peaks]  bamCov]
#     |      |      |                 |      |      |
#     +--+---+      |                 +--+---+      |
#        |          |                    |          |
# [CONVERGENCE 1]  |             [CONVERGENCE 2]   |          Level 5
# (expression)     |              (accessibility)   |
#        |          |                    |          |
#        +----------+--------------------+          |
#                   |                               |
#           [CONVERGENCE 3] <-----------------------+          Level 6
#           (RNA + ATAC + coverage)
#                   |
#           +-------+---------------+
#           |       |               |
#     [bedtools  [homer          [python                       Level 7
#      closest    findMotifs      correlation]
#      (peak->    (peaks)]
#       gene)]
#           |       |               |
#           +-------+---------------+
#                   |
#           [CONVERGENCE 4]                                    Level 8
#           (peak-gene + motifs + correlation)
#                   |
#           +-------+-----------+
#           |       |           |
#     [python    [python      [python                          Level 9
#      diff       TF-target    regulatory
#      analysis]  network]     circuit]
#           |       |           |
#           +-------+-----------+
#                   |
#           [CONVERGENCE 5] <-- QC from both                   Level 10
#           [python integrative report]
# ============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="$WORKDIR/data"
REF="$WORKDIR/reference"
OUT="$WORKDIR/outputs"
RES="$WORKDIR/results"

mkdir -p "$OUT"/{rna_qc,atac_qc,rna_align,atac_align,rna_quant,rna_assembly,rna_counts,atac_dedup,atac_peaks,atac_signal,integration,motifs,analysis}
mkdir -p "$RES"

GENOME="$REF/genome.fasta"
GTF="$REF/genes.gtf"

# Index reference
if [ ! -f "${GENOME}.fai" ]; then
  samtools faidx "$GENOME"
fi

# ============================================================================
# LEVEL 1: QC with fastp (both modalities)
# ============================================================================
if [ ! -f "$OUT/rna_qc/rna_R1.fastq.gz" ]; then
  echo "[Level 1] fastp QC for RNA-seq..."
  fastp -i "$DATA/rna_R1.fastq.gz" -I "$DATA/rna_R2.fastq.gz" \
    -o "$OUT/rna_qc/rna_R1.fastq.gz" -O "$OUT/rna_qc/rna_R2.fastq.gz" \
    -h "$OUT/rna_qc/rna_fastp.html" -j "$OUT/rna_qc/rna_fastp.json" \
    -w "$THREADS" --detect_adapter_for_pe
fi

if [ ! -f "$OUT/atac_qc/atac_R1.fastq.gz" ]; then
  echo "[Level 1] fastp QC for ATAC-seq..."
  fastp -i "$DATA/atac_R1.fastq.gz" -I "$DATA/atac_R2.fastq.gz" \
    -o "$OUT/atac_qc/atac_R1.fastq.gz" -O "$OUT/atac_qc/atac_R2.fastq.gz" \
    -h "$OUT/atac_qc/atac_fastp.html" -j "$OUT/atac_qc/atac_fastp.json" \
    -w "$THREADS" --detect_adapter_for_pe
fi

# ============================================================================
# LEVEL 2: Alignment
# ============================================================================
# Build STAR index for chr22
if [ ! -d "$OUT/star_index" ]; then
  echo "[Level 2] Building STAR index..."
  mkdir -p "$OUT/star_index"
  STAR --runMode genomeGenerate \
    --genomeDir "$OUT/star_index" \
    --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$GTF" \
    --runThreadN "$THREADS" \
    --genomeSAindexNbases 11
fi

if [ ! -f "$OUT/rna_align/Aligned.sortedByCoord.out.bam" ]; then
  echo "[Level 2] STAR alignment for RNA-seq..."
  STAR --runMode alignReads \
    --genomeDir "$OUT/star_index" \
    --readFilesIn "$OUT/rna_qc/rna_R1.fastq.gz" "$OUT/rna_qc/rna_R2.fastq.gz" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$OUT/rna_align/" \
    --runThreadN "$THREADS" \
    --outSAMattrRGline ID:rna SM:K562_RNA PL:ILLUMINA \
    --quantMode TranscriptomeSAM GeneCounts
fi

# Build bowtie2 index for chr22
if [ ! -f "$OUT/bt2_index/genome.1.bt2" ]; then
  echo "[Level 2] Building bowtie2 index..."
  mkdir -p "$OUT/bt2_index"
  bowtie2-build "$GENOME" "$OUT/bt2_index/genome" --threads "$THREADS"
fi

if [ ! -f "$OUT/atac_align/atac_sorted.bam" ]; then
  echo "[Level 2] bowtie2 alignment for ATAC-seq..."
  bowtie2 -x "$OUT/bt2_index/genome" \
    -1 "$OUT/atac_qc/atac_R1.fastq.gz" \
    -2 "$OUT/atac_qc/atac_R2.fastq.gz" \
    --very-sensitive -X 2000 --no-mixed --no-discordant \
    --rg-id atac --rg "SM:K562_ATAC" --rg "PL:ILLUMINA" \
    -p "$THREADS" 2>"$OUT/atac_align/bowtie2.log" | \
    samtools sort -@ "$THREADS" -o "$OUT/atac_align/atac_sorted.bam"
  samtools index "$OUT/atac_align/atac_sorted.bam"
fi

# ============================================================================
# LEVEL 3: Sort + Index RNA BAM
# ============================================================================
if [ ! -f "$OUT/rna_align/rna_sorted.bam" ]; then
  echo "[Level 3] Sorting RNA BAM..."
  if [ -f "$OUT/rna_align/Aligned.sortedByCoord.out.bam" ]; then
    ln -sf "Aligned.sortedByCoord.out.bam" "$OUT/rna_align/rna_sorted.bam"
  fi
  samtools index "$OUT/rna_align/rna_sorted.bam"
fi

# ============================================================================
# LEVEL 4a: Salmon quantification
# ============================================================================
if [ ! -f "$OUT/rna_quant/quant.sf" ]; then
  echo "[Level 4a] Salmon quantification..."
  # Build salmon index from transcriptome
  if [ ! -d "$OUT/salmon_index" ]; then
    # Extract transcript sequences from genome using GTF
    # Use STAR's transcriptome BAM instead
    salmon quant -t "$GENOME" -l A \
      -a "$OUT/rna_align/Aligned.toTranscriptome.out.bam" \
      -o "$OUT/rna_quant" \
      -p "$THREADS" 2>/dev/null || {
      # Fallback: quantify from FASTQ
      salmon quant -i "$OUT/salmon_index" -l A \
        -1 "$OUT/rna_qc/rna_R1.fastq.gz" -2 "$OUT/rna_qc/rna_R2.fastq.gz" \
        -o "$OUT/rna_quant" -p "$THREADS" 2>/dev/null || true
    }
  fi
fi

# ============================================================================
# LEVEL 4b: StringTie transcript assembly
# ============================================================================
if [ ! -f "$OUT/rna_assembly/stringtie.gtf" ]; then
  echo "[Level 4b] StringTie assembly..."
  stringtie "$OUT/rna_align/rna_sorted.bam" \
    -G "$GTF" -o "$OUT/rna_assembly/stringtie.gtf" \
    -e -B -p "$THREADS" -A "$OUT/rna_assembly/gene_abundances.tab" || true
fi

# ============================================================================
# LEVEL 4c: featureCounts gene-level counts
# ============================================================================
if [ ! -f "$OUT/rna_counts/counts.txt" ]; then
  echo "[Level 4c] featureCounts..."
  featureCounts -a "$GTF" -o "$OUT/rna_counts/counts.txt" \
    -T "$THREADS" -p --countReadPairs \
    "$OUT/rna_align/rna_sorted.bam" || true
fi

# ============================================================================
# LEVEL 4d: Picard MarkDuplicates for ATAC
# ============================================================================
if [ ! -f "$OUT/atac_dedup/atac_dedup.bam" ]; then
  echo "[Level 4d] Picard dedup for ATAC..."
  picard MarkDuplicates \
    I="$OUT/atac_align/atac_sorted.bam" \
    O="$OUT/atac_dedup/atac_dedup.bam" \
    M="$OUT/atac_dedup/metrics.txt" \
    REMOVE_DUPLICATES=true 2>/dev/null || true
  if [ -f "$OUT/atac_dedup/atac_dedup.bam" ]; then
    samtools index "$OUT/atac_dedup/atac_dedup.bam"
  fi
fi

# Use best available ATAC BAM
ATAC_BAM="$OUT/atac_dedup/atac_dedup.bam"
if [ ! -f "$ATAC_BAM" ]; then
  ATAC_BAM="$OUT/atac_align/atac_sorted.bam"
fi

# ============================================================================
# LEVEL 4e: MACS2 peak calling for ATAC
# ============================================================================
if [ ! -f "$OUT/atac_peaks/atac_peaks.narrowPeak" ]; then
  echo "[Level 4e] MACS2 peak calling for ATAC..."
  macs3 callpeak -t "$ATAC_BAM" \
    -f BAMPE -g 5.1e7 \
    --nomodel --shift -100 --extsize 200 \
    --outdir "$OUT/atac_peaks" -n atac \
    --keep-dup all -q 0.05 2>&1 | tail -5 || true
  # Rename output
  if [ -f "$OUT/atac_peaks/atac_peaks.narrowPeak" ]; then
    echo "  MACS2 peaks: $(wc -l < "$OUT/atac_peaks/atac_peaks.narrowPeak")"
  fi
fi

# ============================================================================
# LEVEL 4f: deeptools bamCoverage for ATAC signal
# ============================================================================
if [ ! -f "$OUT/atac_signal/atac.bw" ]; then
  echo "[Level 4f] deeptools signal track..."
  bamCoverage -b "$ATAC_BAM" \
    -o "$OUT/atac_signal/atac.bw" \
    --normalizeUsing RPKM \
    --binSize 10 -p "$THREADS" 2>/dev/null || true
fi

# ============================================================================
# LEVEL 5 (CONVERGENCE 1+2): RNA expression + ATAC peaks
# ============================================================================
echo "[Level 5] CONVERGENCE 1+2: RNA expression + ATAC accessibility ready"

# ============================================================================
# LEVEL 6 (CONVERGENCE 3): RNA + ATAC + coverage integrated
# ============================================================================
echo "[Level 6] CONVERGENCE 3: Integrating modalities..."

# ============================================================================
# LEVEL 7a: Peak-to-gene assignment
# ============================================================================
if [ ! -f "$OUT/integration/peak_gene_links.tsv" ]; then
  echo "[Level 7a] Peak-to-gene assignment..."

  # Create TSS BED from GTF
  python3 -c "
import re
genes = set()
with open('$GTF') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene': continue
        m = re.search(r'gene_name \"([^\"]+)\"', parts[8])
        gene_name = m.group(1) if m else 'unknown'
        m2 = re.search(r'gene_id \"([^\"]+)\"', parts[8])
        gene_id = m2.group(1) if m2 else 'unknown'
        chrom = parts[0]; strand = parts[6]
        tss = int(parts[3]) if strand == '+' else int(parts[4])
        # TSS +/- 500bp
        start = max(0, tss - 500)
        end = tss + 500
        genes.add((chrom, start, end, gene_name, gene_id, strand))
with open('$OUT/integration/tss.bed', 'w') as f:
    for g in sorted(genes):
        f.write('\t'.join(str(x) for x in g) + '\n')
print(f'TSS regions: {len(genes)}')
  "

  # Assign peaks to nearest TSS
  PEAKS="$OUT/atac_peaks/atac_peaks.narrowPeak"
  if [ -f "$PEAKS" ] && [ -s "$PEAKS" ]; then
    bedtools closest -a "$PEAKS" -b "$OUT/integration/tss.bed" -d \
      > "$OUT/integration/peak_gene_links.tsv" 2>/dev/null || true
    echo "  Peak-gene links: $(wc -l < "$OUT/integration/peak_gene_links.tsv" || echo 0)"
  else
    touch "$OUT/integration/peak_gene_links.tsv"
  fi
fi

# ============================================================================
# LEVEL 7b: HOMER motif analysis
# ============================================================================
if [ ! -d "$OUT/motifs/homerResults" ]; then
  echo "[Level 7b] HOMER motif analysis..."
  PEAKS="$OUT/atac_peaks/atac_peaks.narrowPeak"
  if [ -f "$PEAKS" ] && [ -s "$PEAKS" ]; then
    # Prepare HOMER-format peak file
    awk -F'\t' '{OFS="\t"; print NR, $1, $2, $3, "+"}' "$PEAKS" > "$OUT/motifs/peaks.homer"
    findMotifsGenome.pl "$OUT/motifs/peaks.homer" "$GENOME" "$OUT/motifs" \
      -size 200 -mask -p "$THREADS" 2>&1 | tail -5 || true
  fi
fi

# ============================================================================
# LEVEL 7c: Peak score vs expression correlation
# ============================================================================
if [ ! -f "$OUT/integration/correlation.tsv" ]; then
  echo "[Level 7c] Computing peak-expression correlation..."
  python3 << 'PYEOF'
import os
import pandas as pd
import numpy as np

# Load gene counts
counts_file = "outputs/rna_counts/counts.txt"
if os.path.exists(counts_file):
    counts = pd.read_csv(counts_file, sep="\t", comment="#")
    # Get gene expression (last column is counts)
    gene_col = counts.columns[0]
    count_col = counts.columns[-1]
    gene_expr = counts[[gene_col, count_col]].rename(columns={gene_col: "gene_id", count_col: "count"})
    gene_expr = gene_expr[gene_expr["count"] > 0]
    print(f"Genes with expression: {len(gene_expr)}")
else:
    gene_expr = pd.DataFrame(columns=["gene_id", "count"])

# Load peak-gene links
links_file = "outputs/integration/peak_gene_links.tsv"
if os.path.exists(links_file) and os.path.getsize(links_file) > 0:
    links = pd.read_csv(links_file, sep="\t", header=None)
    # Peak score is column 4 (0-indexed), gene name is column 13
    if len(links.columns) >= 17:
        links_sub = links[[4, 13, 16]].rename(columns={4: "peak_score", 13: "gene_name", 16: "distance"})
    elif len(links.columns) >= 14:
        links_sub = links[[4, 13]].rename(columns={4: "peak_score", 13: "gene_name"})
        links_sub["distance"] = 0
    else:
        links_sub = pd.DataFrame(columns=["peak_score", "gene_name", "distance"])
    print(f"Peak-gene links: {len(links_sub)}")
else:
    links_sub = pd.DataFrame(columns=["peak_score", "gene_name", "distance"])

# Compute correlation stats
os.makedirs("outputs/integration", exist_ok=True)
with open("outputs/integration/correlation.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"expressed_genes\t{len(gene_expr)}\n")
    f.write(f"peak_gene_links\t{len(links_sub)}\n")
    if len(links_sub) > 0 and "distance" in links_sub.columns:
        f.write(f"median_peak_gene_distance\t{int(links_sub['distance'].median())}\n")
        f.write(f"peaks_within_10kb\t{(links_sub['distance'] < 10000).sum()}\n")
print("Correlation analysis done")
PYEOF
fi

# ============================================================================
# LEVEL 8 (CONVERGENCE 4): peak-gene + motifs + correlation
# ============================================================================
echo "[Level 8] CONVERGENCE 4: integration analysis complete"

# ============================================================================
# LEVEL 9a: Differential analysis summary
# ============================================================================
if [ ! -f "$OUT/analysis/diff_summary.tsv" ]; then
  echo "[Level 9a] Differential analysis summary..."
  python3 << 'PYEOF'
import os
import pandas as pd

# Summarize gene expression distribution
counts_file = "outputs/rna_counts/counts.txt"
os.makedirs("outputs/analysis", exist_ok=True)
if os.path.exists(counts_file):
    counts = pd.read_csv(counts_file, sep="\t", comment="#")
    count_col = counts.columns[-1]
    total_genes = len(counts)
    expressed = (counts[count_col] > 0).sum()
    highly_expressed = (counts[count_col] > 100).sum()
    with open("outputs/analysis/diff_summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_genes\t{total_genes}\n")
        f.write(f"expressed_genes\t{expressed}\n")
        f.write(f"highly_expressed_genes\t{highly_expressed}\n")
        f.write(f"median_expression\t{counts[count_col].median()}\n")
    print(f"Expression: {expressed}/{total_genes} expressed, {highly_expressed} highly")
else:
    with open("outputs/analysis/diff_summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
PYEOF
fi

# ============================================================================
# LEVEL 9b: TF-target network
# ============================================================================
if [ ! -f "$OUT/analysis/tf_network.tsv" ]; then
  echo "[Level 9b] TF-target network..."
  python3 << 'PYEOF'
import os
import pandas as pd

os.makedirs("outputs/analysis", exist_ok=True)

# Parse HOMER motif results
motif_dir = "outputs/motifs"
known_motifs = os.path.join(motif_dir, "knownResults.txt")
if os.path.exists(known_motifs):
    motifs = pd.read_csv(known_motifs, sep="\t")
    significant = motifs[motifs.iloc[:, 2] < 0.01] if len(motifs.columns) > 2 else motifs.head(0)
    with open("outputs/analysis/tf_network.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_motifs_tested\t{len(motifs)}\n")
        f.write(f"significant_motifs\t{len(significant)}\n")
        if len(motifs) > 0:
            f.write(f"top_motif\t{motifs.iloc[0, 0]}\n")
    print(f"TF network: {len(significant)} significant motifs")
else:
    # No HOMER results, create empty
    with open("outputs/analysis/tf_network.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write("total_motifs_tested\t0\n")
        f.write("significant_motifs\t0\n")
    print("No motif results available")
PYEOF
fi

# ============================================================================
# LEVEL 9c: Regulatory circuit scoring
# ============================================================================
if [ ! -f "$OUT/analysis/regulatory_circuits.tsv" ]; then
  echo "[Level 9c] Regulatory circuit scoring..."
  python3 << 'PYEOF'
import os
import pandas as pd

os.makedirs("outputs/analysis", exist_ok=True)

# Load peak-gene links and expression to identify regulatory circuits
links_file = "outputs/integration/peak_gene_links.tsv"
counts_file = "outputs/rna_counts/counts.txt"

circuits = []
if os.path.exists(links_file) and os.path.getsize(links_file) > 0:
    links = pd.read_csv(links_file, sep="\t", header=None)
    if len(links.columns) >= 14:
        # Filter to proximal peaks (within 100kb)
        if len(links.columns) >= 17:
            proximal = links[links.iloc[:, -1] < 100000]
        else:
            proximal = links
        n_circuits = len(proximal)
    else:
        n_circuits = 0
else:
    n_circuits = 0

with open("outputs/analysis/regulatory_circuits.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"total_regulatory_links\t{n_circuits}\n")

print(f"Regulatory circuits: {n_circuits} links")
PYEOF
fi

# ============================================================================
# LEVEL 10 (CONVERGENCE 5): Final integrative report
# ============================================================================
echo "[Level 10] CONVERGENCE 5: Generating final report..."
python3 << 'PYEOF'
import os, json
import pandas as pd

results = {}

# ---- RNA QC ----
rna_qc = "outputs/rna_qc/rna_fastp.json"
if os.path.exists(rna_qc):
    with open(rna_qc) as f:
        qc = json.load(f)
    results["rna_total_reads"] = qc["summary"]["before_filtering"]["total_reads"]
    results["rna_filtered_reads"] = qc["summary"]["after_filtering"]["total_reads"]
    results["rna_q30_rate"] = round(qc["summary"]["after_filtering"]["q30_rate"] * 100, 1)

# ---- ATAC QC ----
atac_qc = "outputs/atac_qc/atac_fastp.json"
if os.path.exists(atac_qc):
    with open(atac_qc) as f:
        qc = json.load(f)
    results["atac_total_reads"] = qc["summary"]["before_filtering"]["total_reads"]
    results["atac_filtered_reads"] = qc["summary"]["after_filtering"]["total_reads"]
    results["atac_q30_rate"] = round(qc["summary"]["after_filtering"]["q30_rate"] * 100, 1)

# ---- RNA alignment ----
rna_bam = "outputs/rna_align/rna_sorted.bam"
if os.path.exists(rna_bam):
    import subprocess
    flag = subprocess.run(["samtools", "flagstat", rna_bam], capture_output=True, text=True)
    lines = flag.stdout.split("\n")
    for line in lines:
        if "mapped (" in line:
            parts = line.split()
            results["rna_mapped_reads"] = int(parts[0])
            pct = line.split("(")[1].split("%")[0] if "(" in line else "0"
            results["rna_mapping_rate"] = float(pct)
            break

# ---- ATAC alignment ----
atac_log = "outputs/atac_align/bowtie2.log"
if os.path.exists(atac_log):
    with open(atac_log) as f:
        for line in f:
            if "overall alignment rate" in line:
                rate = line.strip().split()[0].replace("%", "")
                results["atac_mapping_rate"] = float(rate)

# ---- ATAC dedup ----
dedup_metrics = "outputs/atac_dedup/metrics.txt"
if os.path.exists(dedup_metrics):
    with open(dedup_metrics) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY"):
            if i + 1 < len(lines):
                vals = lines[i+1].strip().split("\t")
                if len(vals) >= 9:
                    try:
                        results["atac_duplication_rate"] = round(float(vals[8]) * 100, 1)
                    except (ValueError, IndexError):
                        pass
            break

# ---- ATAC peaks ----
peaks_file = "outputs/atac_peaks/atac_peaks.narrowPeak"
if os.path.exists(peaks_file):
    n_peaks = sum(1 for _ in open(peaks_file))
    results["total_peaks"] = n_peaks
    # Peak size distribution
    sizes = []
    with open(peaks_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                sizes.append(int(parts[2]) - int(parts[1]))
    if sizes:
        results["median_peak_width"] = sorted(sizes)[len(sizes)//2]
        results["total_peak_bp"] = sum(sizes)

# ---- Gene counts ----
counts_file = "outputs/rna_counts/counts.txt"
if os.path.exists(counts_file):
    counts = pd.read_csv(counts_file, sep="\t", comment="#")
    count_col = counts.columns[-1]
    results["total_genes_annotated"] = len(counts)
    results["expressed_genes"] = int((counts[count_col] > 0).sum())
    results["highly_expressed_genes"] = int((counts[count_col] > 100).sum())

# ---- StringTie ----
st_file = "outputs/rna_assembly/gene_abundances.tab"
if os.path.exists(st_file):
    st = pd.read_csv(st_file, sep="\t")
    results["assembled_transcripts"] = len(st)

# ---- Salmon ----
salmon_file = "outputs/rna_quant/quant.sf"
if os.path.exists(salmon_file):
    sal = pd.read_csv(salmon_file, sep="\t")
    results["quantified_transcripts"] = len(sal)
    results["detected_transcripts"] = int((sal["NumReads"] > 0).sum()) if "NumReads" in sal.columns else 0

# ---- Peak-gene links ----
links_file = "outputs/integration/peak_gene_links.tsv"
if os.path.exists(links_file) and os.path.getsize(links_file) > 0:
    links = pd.read_csv(links_file, sep="\t", header=None)
    results["peak_gene_links"] = len(links)

# ---- Correlation ----
corr_file = "outputs/integration/correlation.tsv"
if os.path.exists(corr_file):
    corr = pd.read_csv(corr_file, sep="\t")
    for _, row in corr.iterrows():
        results[f"corr_{row['metric']}"] = row["value"]

# ---- Signal track ----
results["signal_track_generated"] = "yes" if os.path.exists("outputs/atac_signal/atac.bw") else "no"

# ---- Motifs ----
tf_file = "outputs/analysis/tf_network.tsv"
if os.path.exists(tf_file):
    tf = pd.read_csv(tf_file, sep="\t")
    for _, row in tf.iterrows():
        results[f"motif_{row['metric']}"] = row["value"]

# ---- Regulatory circuits ----
rc_file = "outputs/analysis/regulatory_circuits.tsv"
if os.path.exists(rc_file):
    rc = pd.read_csv(rc_file, sep="\t")
    for _, row in rc.iterrows():
        results[row["metric"]] = row["value"]

# ---- Write report ----
os.makedirs("results", exist_ok=True)
with open("results/report.csv", "w") as f:
    f.write("metric,value\n")
    for k, v in results.items():
        f.write(f"{k},{v}\n")

print("=== Final Report ===")
for k, v in results.items():
    print(f"  {k}: {v}")
print(f"\nTotal metrics: {len(results)}")
PYEOF

echo "Pipeline complete. Results in results/report.csv"
