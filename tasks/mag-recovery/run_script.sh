#!/usr/bin/env bash
set -euo pipefail

# MAG Recovery from Metagenomic Reads
# Data: Human gut metagenome (DRR256963, Illumina HiSeq 2500, PE 150bp)
#
# DAG Structure (depth=10, convergence=4):
#
#   [R1.fq + R2.fq]
#       │
#     fastp(QC)                                    (Step 1)
#       │
#     MEGAHIT(assembly)                            (Step 2)
#       │
#   ┌───┴─────────────────────────────┐
#   │                                 │
# quast(stats)                  bowtie2-build      (Step 3: parallel)
#   │                                 │
#   │                           bowtie2 align      (Step 4)
#   │                                 │
#   │                         samtools sort+idx    (Step 5)
#   │                                 │
#   │                         depth calculation    (Step 6)
#   │                                 │
#   │                   ┌─────────────┼──────────┐
#   │               MetaBAT2      MaxBin2     prodigal  (Step 7: parallel)
#   │                   └─────────────┼──────────┘
#   │                                 │
#   │                           DAS_Tool            (Step 8: CONVERGE #1)
#   │                                 │
#   └─────────────────────────────────┤
#                                     │
#                              bin assessment       (Step 9: CONVERGE #2 - assembly+bins)
#                          ┌──────────┼──────────┐
#                       per-bin    gene counts   GC/size  (parallel)
#                       N50/stats  per bin       analysis
#                          └──────────┼──────────┘
#                                     │
#                                final report       (Step 10: CONVERGE #3+#4)

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="${WORKDIR}/data"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}/qc" "${OUT}/assembly" "${OUT}/mapping" "${OUT}/binning"
mkdir -p "${OUT}/metabat2" "${OUT}/maxbin2" "${OUT}/dastool" "${OUT}/quast"
mkdir -p "${OUT}/gene_prediction" "${RESULTS}"

R1="${DATA}/reads_R1.fastq.gz"
R2="${DATA}/reads_R2.fastq.gz"

# ============================================================================
# Step 1: Read QC with fastp
# ============================================================================
echo "[Step 1] Running fastp QC..."
if [ ! -f "${OUT}/qc/trimmed_R1.fastq.gz" ]; then
    fastp \
        -i "$R1" -I "$R2" \
        -o "${OUT}/qc/trimmed_R1.fastq.gz" \
        -O "${OUT}/qc/trimmed_R2.fastq.gz" \
        --json "${OUT}/qc/fastp.json" \
        --html "${OUT}/qc/fastp.html" \
        --thread "$THREADS" \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --detect_adapter_for_pe \
        2>&1 | tail -5
fi
echo "[Step 1] QC done."

# ============================================================================
# Step 2: Metagenomic assembly with MEGAHIT
# ============================================================================
echo "[Step 2] Running MEGAHIT assembly..."
if [ ! -f "${OUT}/assembly/final.contigs.fa" ]; then
    rm -rf "${OUT}/assembly"
    megahit \
        -1 "${OUT}/qc/trimmed_R1.fastq.gz" \
        -2 "${OUT}/qc/trimmed_R2.fastq.gz" \
        -o "${OUT}/assembly" \
        -t "$THREADS" \
        --min-contig-len 1000 \
        2>&1 | tail -10
fi
echo "[Step 2] Assembly done."

# Filter to contigs >= 1500bp for binning
echo "  Filtering contigs >= 1500bp..."
python3 << 'FILTER_PY'
import os
out = os.environ.get("OUT", "outputs")
contigs = f"{out}/assembly/final.contigs.fa"
filtered = f"{out}/assembly/contigs_1500.fa"
with open(contigs) as f, open(filtered, 'w') as fo:
    seq = ""
    header = ""
    for line in f:
        if line.startswith(">"):
            if header and len(seq) >= 1500:
                fo.write(header + seq + "\n")
            header = line
            seq = ""
        else:
            seq += line.strip()
    if header and len(seq) >= 1500:
        fo.write(header + seq + "\n")
FILTER_PY
CONTIGS="${OUT}/assembly/contigs_1500.fa"
N_CONTIGS=$(grep -c "^>" "$CONTIGS" || true)
echo "  Contigs >= 1500bp: ${N_CONTIGS}"

# ============================================================================
# Step 3a: Assembly quality with QUAST (parallel with 3b)
# ============================================================================
echo "[Step 3a] Running QUAST..."
if [ ! -f "${OUT}/quast/report.tsv" ]; then
    quast "$CONTIGS" \
        -o "${OUT}/quast" \
        --threads "$THREADS" \
        --min-contig 1000 \
        2>&1 | tail -5
fi

# ============================================================================
# Step 3b: Build bowtie2 index
# ============================================================================
echo "[Step 3b] Building bowtie2 index..."
if [ ! -f "${OUT}/mapping/contigs.1.bt2" ]; then
    bowtie2-build \
        --threads "$THREADS" \
        "$CONTIGS" \
        "${OUT}/mapping/contigs" \
        2>&1 | tail -3
fi

# ============================================================================
# Step 4: Map reads back to contigs
# ============================================================================
echo "[Step 4] Mapping reads to contigs..."
if [ ! -f "${OUT}/mapping/mapped.bam" ]; then
    bowtie2 \
        -x "${OUT}/mapping/contigs" \
        -1 "${OUT}/qc/trimmed_R1.fastq.gz" \
        -2 "${OUT}/qc/trimmed_R2.fastq.gz" \
        --threads "$THREADS" \
        --very-sensitive \
        -S "${OUT}/mapping/mapped.sam" \
        2>&1 | tail -5
    samtools view -bS -@ "$THREADS" "${OUT}/mapping/mapped.sam" > "${OUT}/mapping/mapped.bam"
    rm -f "${OUT}/mapping/mapped.sam"
fi

# ============================================================================
# Step 5: Sort and index BAM
# ============================================================================
echo "[Step 5] Sorting and indexing BAM..."
if [ ! -f "${OUT}/mapping/sorted.bam.bai" ]; then
    samtools sort -@ "$THREADS" -o "${OUT}/mapping/sorted.bam" "${OUT}/mapping/mapped.bam"
    samtools index "${OUT}/mapping/sorted.bam"
    rm -f "${OUT}/mapping/mapped.bam"
fi
echo "  Mapping stats:"
samtools flagstat "${OUT}/mapping/sorted.bam" | head -5

# ============================================================================
# Step 6: Calculate depth
# ============================================================================
echo "[Step 6] Calculating contig depth..."
if [ ! -f "${OUT}/mapping/depth.txt" ]; then
    jgi_summarize_bam_contig_depths \
        --outputDepth "${OUT}/mapping/depth.txt" \
        "${OUT}/mapping/sorted.bam" \
        2>&1 | tail -3
fi

# ============================================================================
# Step 7a: MetaBAT2 binning (parallel)
# ============================================================================
echo "[Step 7a] Running MetaBAT2..."
if [ ! -d "${OUT}/metabat2/bins" ] || [ -z "$(ls ${OUT}/metabat2/bins/ 2>/dev/null)" ]; then
    mkdir -p "${OUT}/metabat2/bins"
    metabat2 \
        -i "$CONTIGS" \
        -a "${OUT}/mapping/depth.txt" \
        -o "${OUT}/metabat2/bins/bin" \
        --minContig 1500 \
        -t "$THREADS" \
        2>&1 | tail -3
fi
METABAT_BINS=$(ls "${OUT}/metabat2/bins/"*.fa 2>/dev/null | wc -l || echo 0)
echo "  MetaBAT2 bins: ${METABAT_BINS}"

# ============================================================================
# Step 7b: MaxBin2 binning (parallel)
# ============================================================================
echo "[Step 7b] Running MaxBin2..."
if [ ! -d "${OUT}/maxbin2/bins" ] || [ -z "$(ls ${OUT}/maxbin2/bins/ 2>/dev/null)" ]; then
    mkdir -p "${OUT}/maxbin2/bins"
    # Create abundance file from depth
    tail -n +2 "${OUT}/mapping/depth.txt" | awk -F'\t' '{print $1"\t"$3}' > "${OUT}/maxbin2/abundance.txt"
    run_MaxBin.pl \
        -contig "$CONTIGS" \
        -abund "${OUT}/maxbin2/abundance.txt" \
        -out "${OUT}/maxbin2/bins/bin" \
        -thread "$THREADS" \
        2>&1 | tail -5 || true
fi
MAXBIN_BINS=$(ls "${OUT}/maxbin2/bins/"*.fasta 2>/dev/null | wc -l || echo 0)
echo "  MaxBin2 bins: ${MAXBIN_BINS}"

# ============================================================================
# Step 7c: Prodigal gene prediction on contigs (parallel, needed by DAS_Tool)
# ============================================================================
echo "[Step 7c] Running Prodigal..."
if [ ! -f "${OUT}/gene_prediction/proteins.faa" ]; then
    prodigal \
        -i "$CONTIGS" \
        -o "${OUT}/gene_prediction/genes.gff" \
        -a "${OUT}/gene_prediction/proteins.faa" \
        -p meta \
        -f gff \
        2>&1 | tail -3
fi
N_GENES=$(grep -c "^>" "${OUT}/gene_prediction/proteins.faa" || echo 0)
echo "  Predicted genes: ${N_GENES}"

# ============================================================================
# Step 8: DAS_Tool bin refinement — CONVERGE #1 (multiple binners)
# ============================================================================
echo "[Step 8] Running DAS_Tool (convergence #1)..."

# Create scaffold-to-bin files for DAS_Tool
echo "  Creating scaffold-to-bin mappings..."
# MetaBAT2
> "${OUT}/dastool/metabat2_s2b.tsv"
for BIN in "${OUT}/metabat2/bins/"*.fa; do
    [ -f "$BIN" ] || continue
    BNAME=$(basename "$BIN" .fa)
    grep "^>" "$BIN" | sed "s/>//" | awk -v b="$BNAME" '{print $1"\t"b}' >> "${OUT}/dastool/metabat2_s2b.tsv"
done

# MaxBin2
> "${OUT}/dastool/maxbin2_s2b.tsv"
for BIN in "${OUT}/maxbin2/bins/"*.fasta; do
    [ -f "$BIN" ] || continue
    BNAME=$(basename "$BIN" .fasta)
    grep "^>" "$BIN" | sed "s/>//" | awk -v b="$BNAME" '{print $1"\t"b}' >> "${OUT}/dastool/maxbin2_s2b.tsv"
done

# Run DAS_Tool
DASTOOL_INPUT=""
DASTOOL_LABELS=""
if [ -s "${OUT}/dastool/metabat2_s2b.tsv" ]; then
    DASTOOL_INPUT="${OUT}/dastool/metabat2_s2b.tsv"
    DASTOOL_LABELS="metabat2"
fi
if [ -s "${OUT}/dastool/maxbin2_s2b.tsv" ]; then
    if [ -n "$DASTOOL_INPUT" ]; then
        DASTOOL_INPUT="${DASTOOL_INPUT},${OUT}/dastool/maxbin2_s2b.tsv"
        DASTOOL_LABELS="${DASTOOL_LABELS},maxbin2"
    else
        DASTOOL_INPUT="${OUT}/dastool/maxbin2_s2b.tsv"
        DASTOOL_LABELS="maxbin2"
    fi
fi

if [ -n "$DASTOOL_INPUT" ] && [ ! -f "${OUT}/dastool/results_DASTool_summary.tsv" ]; then
    DAS_Tool \
        -i "$DASTOOL_INPUT" \
        -l "$DASTOOL_LABELS" \
        -c "$CONTIGS" \
        -o "${OUT}/dastool/results" \
        --threads "$THREADS" \
        --search_engine diamond \
        --write_bins \
        --proteins "${OUT}/gene_prediction/proteins.faa" \
        2>&1 | tail -10 || true
fi

# Count final bins
FINAL_BINS=0
if [ -d "${OUT}/dastool/results_DASTool_bins" ]; then
    FINAL_BINS=$(ls "${OUT}/dastool/results_DASTool_bins/"*.fa 2>/dev/null | wc -l || echo 0)
fi
# If DAS_Tool didn't produce bins, fall back to metabat2
if [ "$FINAL_BINS" -eq 0 ]; then
    echo "  DAS_Tool produced no bins, using MetaBAT2 bins directly"
    mkdir -p "${OUT}/dastool/results_DASTool_bins"
    if [ "$METABAT_BINS" -gt 0 ]; then
        cp "${OUT}/metabat2/bins/"*.fa "${OUT}/dastool/results_DASTool_bins/" 2>/dev/null || true
        FINAL_BINS=$METABAT_BINS
    fi
fi
echo "  Final refined bins: ${FINAL_BINS}"

# ============================================================================
# Step 9: Bin quality assessment — CONVERGE #2 (assembly stats + bin stats)
# ============================================================================
echo "[Step 9] Assessing bin quality (convergence #2)..."
python3 << 'BIN_QC_PY'
import os, glob

out = os.environ.get("OUT", "outputs")
bin_dir = f"{out}/dastool/results_DASTool_bins"

bins = glob.glob(f"{bin_dir}/*.fa")
bin_stats = []

for bpath in sorted(bins):
    bname = os.path.basename(bpath).replace('.fa', '')
    seqs = []
    with open(bpath) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    seqs.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq:
            seqs.append(seq)

    total_len = sum(len(s) for s in seqs)
    n_contigs = len(seqs)
    gc = sum(s.count('G') + s.count('C') + s.count('g') + s.count('c') for s in seqs) / max(total_len, 1) * 100

    # N50
    lengths = sorted([len(s) for s in seqs], reverse=True)
    cumsum = 0
    n50 = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total_len / 2:
            n50 = l
            break

    bin_stats.append({
        'name': bname,
        'total_length': total_len,
        'n_contigs': n_contigs,
        'n50': n50,
        'gc_pct': round(gc, 2),
        'largest_contig': max(lengths) if lengths else 0,
    })

# Write bin stats
with open(f"{out}/dastool/bin_quality.tsv", 'w') as f:
    f.write("bin\ttotal_length\tn_contigs\tn50\tgc_pct\tlargest_contig\n")
    for b in bin_stats:
        f.write(f"{b['name']}\t{b['total_length']}\t{b['n_contigs']}\t{b['n50']}\t{b['gc_pct']}\t{b['largest_contig']}\n")
        print(f"  {b['name']}: {b['total_length']}bp, {b['n_contigs']} contigs, N50={b['n50']}, GC={b['gc_pct']}%")

print(f"Total bins: {len(bin_stats)}")
BIN_QC_PY

# ============================================================================
# Step 10: Generate final report — CONVERGE #3+#4 (assembly + binning + genes)
# ============================================================================
echo "[Step 10] Generating final report..."
python3 << 'REPORT_PY'
import json, os, glob

out = os.environ.get("OUT", "outputs")
results = os.environ.get("RESULTS", "results")

report = []

# fastp stats
fastp_json = f"{out}/qc/fastp.json"
if os.path.exists(fastp_json):
    with open(fastp_json) as f:
        fj = json.load(f)
    report.append(("total_reads_before", fj["summary"]["before_filtering"]["total_reads"]))
    report.append(("total_reads_after", fj["summary"]["after_filtering"]["total_reads"]))
    report.append(("q30_rate_before", round(fj["summary"]["before_filtering"]["q30_rate"] * 100, 2)))
    report.append(("q30_rate_after", round(fj["summary"]["after_filtering"]["q30_rate"] * 100, 2)))

# QUAST stats
quast_report = f"{out}/quast/report.tsv"
if os.path.exists(quast_report):
    with open(quast_report) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                metric = parts[0].strip()
                value = parts[1].strip()
                if metric == "# contigs (>= 0 bp)":
                    report.append(("total_contigs", value))
                elif metric == "Total length (>= 0 bp)":
                    report.append(("total_assembly_length", value))
                elif metric == "N50":
                    report.append(("assembly_n50", value))
                elif metric == "Largest contig":
                    report.append(("largest_contig", value))
                elif metric == "GC (%)":
                    report.append(("assembly_gc_pct", value))

# Mapping stats
flagstat = os.popen(f"samtools flagstat {out}/mapping/sorted.bam 2>/dev/null").read()
for line in flagstat.split('\n'):
    if 'mapped' in line and 'primary' not in line and '%' in line:
        parts = line.split()
        report.append(("mapped_reads", parts[0]))
        pct = line.split('(')[1].split('%')[0] if '(' in line else '0'
        report.append(("mapping_pct", pct.strip()))
        break

# Gene prediction
genes_file = f"{out}/gene_prediction/proteins.faa"
if os.path.exists(genes_file):
    n_genes = sum(1 for l in open(genes_file) if l.startswith(">"))
    report.append(("predicted_genes", n_genes))

# Binning results
metabat_bins = len(glob.glob(f"{out}/metabat2/bins/*.fa"))
maxbin_bins = len(glob.glob(f"{out}/maxbin2/bins/*.fasta"))
final_bins = len(glob.glob(f"{out}/dastool/results_DASTool_bins/*.fa"))

report.append(("metabat2_bins", metabat_bins))
report.append(("maxbin2_bins", maxbin_bins))
report.append(("refined_bins", final_bins))

# Per-bin stats
bin_quality = f"{out}/dastool/bin_quality.tsv"
if os.path.exists(bin_quality):
    with open(bin_quality) as f:
        header = f.readline()
        bins = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                bins.append({
                    'name': parts[0],
                    'length': int(parts[1]),
                    'contigs': int(parts[2]),
                    'n50': int(parts[3]),
                    'gc': float(parts[4]),
                })

    if bins:
        # Largest bin
        largest = max(bins, key=lambda b: b['length'])
        report.append(("largest_bin_name", largest['name']))
        report.append(("largest_bin_length", largest['length']))
        report.append(("largest_bin_contigs", largest['contigs']))
        report.append(("largest_bin_n50", largest['n50']))
        report.append(("largest_bin_gc_pct", largest['gc']))

        # Summary stats
        total_binned = sum(b['length'] for b in bins)
        report.append(("total_binned_length", total_binned))
        report.append(("mean_bin_size", round(total_binned / len(bins))))
        report.append(("mean_bin_gc_pct", round(sum(b['gc'] for b in bins) / len(bins), 2)))
    else:
        report.append(("largest_bin_name", "none"))
        report.append(("largest_bin_length", 0))
        report.append(("largest_bin_contigs", 0))
        report.append(("largest_bin_n50", 0))
        report.append(("largest_bin_gc_pct", 0))
        report.append(("total_binned_length", 0))
        report.append(("mean_bin_size", 0))
        report.append(("mean_bin_gc_pct", 0))

# DAS_Tool summary
dastool_summary = f"{out}/dastool/results_DASTool_summary.tsv"
if os.path.exists(dastool_summary):
    with open(dastool_summary) as f:
        for line in f:
            if line.startswith("bin\t"):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                report.append(("dastool_score_" + parts[0], parts[-1]))
                break

with open(f"{results}/report.csv", 'w') as f:
    f.write("metric,value\n")
    for m, v in report:
        f.write(f"{m},{v}\n")

print(f"Report: {len(report)} metrics")
for m, v in report:
    print(f"  {m}: {v}")
REPORT_PY

echo ""
echo "========================================="
echo "  MAG Recovery Complete!"
echo "========================================="
echo ""
cat "${RESULTS}/report.csv"
