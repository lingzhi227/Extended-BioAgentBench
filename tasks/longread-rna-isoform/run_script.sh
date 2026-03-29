#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Long-Read RNA Isoform Discovery — DAG (depth=10, convergence=4)
# ============================================================
#
#  reads.fastq.gz (Nanopore direct RNA)
#           │
#     [NanoPlot QC] ────────────────────────── Level 1
#           │
#     [minimap2 -ax splice] ◄── chr22.fa       Level 2
#           │
#     [samtools sort + index]                   Level 3
#           │
#     ┌─────┼──────────────────┐
#     │     │                  │
#  [StringTie              [samtools            Level 4
#   -L long-read]           flagstat]
#     │                        │
#     └────────┬───────────────┘
#              │
#      [CONVERGENCE 1]                          Level 5
#      [gffcompare vs reference GTF]
#              │
#      ┌───────┼──────────────┐
#      │       │              │
#  [python  [bedtools       [python             Level 6
#   isoform  gene body       transcript
#   classify] coverage]      stats]
#      │       │              │
#      └───────┼──────────────┘
#              │
#      [CONVERGENCE 2]                          Level 7
#      (classification + coverage + stats)
#              │
#      ┌───────┼──────────────┐
#      │       │              │
#  [salmon   [python         [python            Level 8
#   quant     expression      novel gene
#   --ont]    filter]         detection]
#      │       │              │
#      └───────┼──────────────┘
#              │
#      [CONVERGENCE 3]                          Level 9
#      (quant + expression + novel)
#              │
#      [CONVERGENCE 4] ◄── NanoPlot QC          Level 10
#      [python report]
#
# Convergence points:
#   C1: StringTie assembly + alignment stats → gffcompare
#   C2: isoform classification + coverage + transcript stats
#   C3: quantification + expression + novel gene detection
#   C4: final report with QC + all analysis results
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORK=$(pwd)
DATA="${WORK}/data"
REF="${WORK}/reference"
OUT="${WORK}/outputs"
RESULTS="${WORK}/results"

mkdir -p "${OUT}"/{qc,alignment,assembly,comparison,coverage,quantification,analysis} "${RESULTS}"

# ─── Level 1: NanoPlot read QC ───
if [ ! -f "${OUT}/qc/NanoStats.txt" ]; then
  echo "[Level 1] Running NanoPlot QC..."
  NanoPlot \
    --fastq "${DATA}/reads.fastq.gz" \
    --outdir "${OUT}/qc" \
    --threads ${THREADS} \
    --loglength \
    --no_static
fi

# Extract QC stats from NanoPlot
TOTAL_READS=$(grep "Number of reads:" "${OUT}/qc/NanoStats.txt" | awk '{gsub(",","",$NF); print $NF}')
MEAN_LENGTH=$(grep "Mean read length:" "${OUT}/qc/NanoStats.txt" | awk '{gsub(",","",$NF); print $NF}')
MEAN_QUAL=$(grep "Mean read quality:" "${OUT}/qc/NanoStats.txt" | awk '{print $NF}')
MEDIAN_LENGTH=$(grep "Median read length:" "${OUT}/qc/NanoStats.txt" | awk '{gsub(",","",$NF); print $NF}')
N50_READ=$(grep "Read length N50:" "${OUT}/qc/NanoStats.txt" | awk '{gsub(",","",$NF); print $NF}')
TOTAL_BASES=$(grep "Total bases:" "${OUT}/qc/NanoStats.txt" | awk '{gsub(",","",$NF); print $NF}')
echo "  Total reads: ${TOTAL_READS}, mean length: ${MEAN_LENGTH}, mean qual: ${MEAN_QUAL}"

# ─── Level 2: minimap2 splice-aware alignment ───
if [ ! -f "${OUT}/alignment/aligned.bam" ]; then
  echo "[Level 2] Running minimap2 splice-aware alignment..."
  minimap2 \
    -ax splice \
    -k14 \
    --secondary=no \
    -t ${THREADS} \
    "${REF}/chr22.fa" \
    "${DATA}/reads.fastq.gz" \
    2> "${OUT}/alignment/minimap2.log" \
  | samtools view -bS -F 4 - \
  > "${OUT}/alignment/aligned.bam"
fi

# ─── Level 3: samtools sort + index ───
if [ ! -f "${OUT}/alignment/aligned.sorted.bam.bai" ]; then
  echo "[Level 3] Sorting and indexing BAM..."
  samtools sort -@ ${THREADS} -o "${OUT}/alignment/aligned.sorted.bam" "${OUT}/alignment/aligned.bam"
  samtools index "${OUT}/alignment/aligned.sorted.bam"
fi

# ─── Level 4a: StringTie long-read transcript assembly ───
if [ ! -f "${OUT}/assembly/stringtie.gtf" ]; then
  echo "[Level 4a] Running StringTie long-read assembly..."
  stringtie \
    "${OUT}/alignment/aligned.sorted.bam" \
    -o "${OUT}/assembly/stringtie.gtf" \
    -G "${REF}/chr22.gtf" \
    -L \
    -p ${THREADS} \
    -v 2> "${OUT}/assembly/stringtie.log"
fi

# ─── Level 4b: samtools flagstat (parallel with assembly) ───
samtools flagstat "${OUT}/alignment/aligned.sorted.bam" > "${OUT}/alignment/flagstat.txt"
MAPPED_READS=$(grep "primary mapped" "${OUT}/alignment/flagstat.txt" | awk '{print $1}')
# Mapping rate relative to total input reads (BAM only has mapped reads due to -F 4)
MAPPING_PCT=$(python3 -c "print(round(${MAPPED_READS}/${TOTAL_READS}*100 if float('${TOTAL_READS}')>0 else 0, 2))")
echo "  Mapped to chr22: ${MAPPED_READS}/${TOTAL_READS} (${MAPPING_PCT}%)"

# ─── Level 5: CONVERGENCE 1 — gffcompare vs reference annotation ───
if [ ! -f "${OUT}/comparison/gffcmp.stats" ]; then
  echo "[Level 5 / CONVERGENCE 1] Running gffcompare..."
  gffcompare \
    -r "${REF}/chr22.gtf" \
    -o "${OUT}/comparison/gffcmp" \
    "${OUT}/assembly/stringtie.gtf"
fi

# Parse gffcompare stats
SENSITIVITY_BASE=$(grep "Base level" "${OUT}/comparison/gffcmp.stats" | head -1 | awk -F'|' '{gsub(/[[:space:]]/,"",$1); sub(/.*:/,"",$1); print $1}')
PRECISION_BASE=$(grep "Base level" "${OUT}/comparison/gffcmp.stats" | head -1 | awk -F'|' '{gsub(/[[:space:]]/,"",$2); print $2}')
TOTAL_TRANSCRIPTS=$(grep "Query mRNAs" "${OUT}/comparison/gffcmp.stats" | awk -F: '{print $2}' | awk '{print $1}')
TOTAL_LOCI=$(grep "Query mRNAs" "${OUT}/comparison/gffcmp.stats" | awk -F'in' '{print $2}' | awk '{print $1}')
echo "  Transcripts: ${TOTAL_TRANSCRIPTS}, Loci: ${TOTAL_LOCI}"
echo "  Base sensitivity: ${SENSITIVITY_BASE}%, precision: ${PRECISION_BASE}%"

# ─── Level 6: Three-way branch — Classification + Coverage + Stats ───
echo "[Level 6] Running isoform classification, coverage, and stats..."

# 6a: Isoform classification
python3 << 'PYEOF'
import os, collections

out = os.environ.get("OUT", "outputs")

# Parse gffcompare annotated GTF for class codes
class_codes = collections.Counter()
transcript_lens = []
exon_counts = collections.Counter()
current_tx = None
current_exons = 0

with open(f"{out}/comparison/gffcmp.annotated.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        attrs = parts[8]

        if feature == "transcript":
            if current_tx:
                exon_counts[current_exons] += 1
            current_exons = 0
            # Extract class_code
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("class_code"):
                    code = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                    class_codes[code] += 1
                    break
            # Extract length
            start = int(parts[3])
            end = int(parts[4])
            transcript_lens.append(end - start)
            current_tx = True
        elif feature == "exon":
            current_exons += 1

    if current_tx:
        exon_counts[current_exons] += 1

# Classification categories
# = : exact match (FSM)
# c : contained (ISM)
# k : containment of reference
# j : novel junction (NIC)
# i : intronic
# u : intergenic (NNC)
# p : polymerase run-on
# x : antisense
# s : sense intronic
# o : other overlap

os.makedirs(f"{out}/analysis", exist_ok=True)
with open(f"{out}/analysis/isoform_classification.tsv", "w") as f:
    f.write("class_code\tdescription\tcount\n")
    desc_map = {"=": "exact_match", "c": "contained", "k": "containment_of_ref",
                "j": "novel_junction", "i": "intronic", "u": "intergenic",
                "p": "polymerase_runon", "x": "antisense", "s": "sense_intronic",
                "o": "other_overlap", "e": "single_exon_partial", ".": "unclassified"}
    for code, count in class_codes.most_common():
        desc = desc_map.get(code, f"other_{code}")
        f.write(f"{code}\t{desc}\t{count}\n")

# Novel isoforms = j + i + u + x + s + o + e
known = class_codes.get("=", 0) + class_codes.get("c", 0)
novel = sum(v for k, v in class_codes.items() if k not in ("=", "c"))
multi_exon = sum(v for k, v in exon_counts.items() if k > 1)

with open(f"{out}/analysis/classification_summary.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"known_isoforms\t{known}\n")
    f.write(f"novel_isoforms\t{novel}\n")
    f.write(f"multi_exon_transcripts\t{multi_exon}\n")
    f.write(f"median_transcript_length\t{sorted(transcript_lens)[len(transcript_lens)//2] if transcript_lens else 0}\n")

print(f"  Classification: {known} known, {novel} novel, {multi_exon} multi-exon")
PYEOF

# 6b: Gene body coverage with bedtools
if [ ! -f "${OUT}/coverage/gene_coverage.tsv" ]; then
  bedtools genomecov \
    -ibam "${OUT}/alignment/aligned.sorted.bam" \
    -d > "${OUT}/coverage/per_base_coverage.tsv"

  # Summarize coverage per gene
  bedtools coverage \
    -a "${REF}/chr22.gtf" \
    -b "${OUT}/alignment/aligned.sorted.bam" \
    -counts \
    | awk -F'\t' '$3=="gene"' \
    | sort -t$'\t' -k10,10nr \
    > "${OUT}/coverage/gene_coverage.tsv"
fi

# 6c: Transcript stats
python3 << 'PYEOF'
import os

out = os.environ.get("OUT", "outputs")

# Parse StringTie GTF for transcript-level stats
transcripts = []
genes = set()
with open(f"{out}/assembly/stringtie.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] == "transcript":
            start = int(parts[3])
            end = int(parts[4])
            length = end - start
            # Extract gene_id and coverage
            attrs = parts[8]
            gene_id = ""
            cov = 0
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1] if '"' in attr else ""
                elif attr.startswith("cov"):
                    try:
                        cov = float(attr.split('"')[1] if '"' in attr else attr.split(" ")[1])
                    except:
                        pass
            genes.add(gene_id)
            transcripts.append({"gene": gene_id, "length": length, "coverage": cov})

# Transcript length distribution
lens = [t["length"] for t in transcripts]
lens.sort()
n50_idx = 0
total_len = sum(lens)
cum = 0
for l in sorted(lens, reverse=True):
    cum += l
    if cum >= total_len / 2:
        n50_idx = l
        break

with open(f"{out}/analysis/transcript_stats.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"total_transcripts\t{len(transcripts)}\n")
    f.write(f"total_genes\t{len(genes)}\n")
    f.write(f"mean_transcript_length\t{round(sum(lens)/len(lens)) if lens else 0}\n")
    f.write(f"median_transcript_length\t{lens[len(lens)//2] if lens else 0}\n")
    f.write(f"transcript_n50\t{n50_idx}\n")
    f.write(f"mean_coverage\t{round(sum(t['coverage'] for t in transcripts)/len(transcripts),2) if transcripts else 0}\n")

print(f"  Transcript stats: {len(transcripts)} transcripts, {len(genes)} genes, N50={n50_idx}")
PYEOF

# ─── Level 7: CONVERGENCE 2 ───
echo "[Level 7 / CONVERGENCE 2] Merging classification, coverage, and transcript stats..."

# ─── Level 8: Three-way branch — Quantification + Expression + Novel detection ───
echo "[Level 8] Running quantification and expression analysis..."

# 8a: Transcript quantification using bedtools coverage
if [ ! -f "${OUT}/quantification/transcript_counts.tsv" ]; then
  # Create BED from StringTie transcripts
  python3 << 'PYEOF'
import os
out = os.environ.get("OUT", "outputs")
os.makedirs(f"{out}/quantification", exist_ok=True)

tx_regions = []
with open(f"{out}/assembly/stringtie.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "transcript":
            continue
        chrom = parts[0]
        start = int(parts[3]) - 1
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        tx_id = ""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("transcript_id"):
                tx_id = attr.split('"')[1] if '"' in attr else ""
                break
        if tx_id:
            tx_regions.append(f"{chrom}\t{start}\t{end}\t{tx_id}\t0\t{strand}\n")

with open(f"{out}/quantification/transcripts.bed", "w") as f:
    f.writelines(tx_regions)
print(f"  Wrote {len(tx_regions)} transcript regions for quantification")
PYEOF

  # Sort BED before coverage
  sort -k1,1 -k2,2n "${OUT}/quantification/transcripts.bed" > "${OUT}/quantification/transcripts.sorted.bed"

  bedtools coverage \
    -a "${OUT}/quantification/transcripts.sorted.bed" \
    -b "${OUT}/alignment/aligned.sorted.bam" \
    -counts -sorted \
    > "${OUT}/quantification/transcript_counts.tsv"
fi

# 8b: Expression analysis
python3 << 'PYEOF'
import os

out = os.environ.get("OUT", "outputs")
os.makedirs(f"{out}/analysis", exist_ok=True)

# Parse bedtools coverage counts
quant_data = []
expressed = 0
total_reads = 0
with open(f"{out}/quantification/transcript_counts.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 7:
            tx_id = parts[3]
            length = int(parts[2]) - int(parts[1])
            count = int(parts[6])
            quant_data.append({"name": tx_id, "count": count, "length": length})
            total_reads += count
            if count > 0:
                expressed += 1

# Compute TPM
for q in quant_data:
    rpk = q["count"] / (q["length"] / 1000) if q["length"] > 0 else 0
    q["rpk"] = rpk
rpk_sum = sum(q["rpk"] for q in quant_data)
for q in quant_data:
    q["tpm"] = round(q["rpk"] / rpk_sum * 1e6, 2) if rpk_sum > 0 else 0

quant_data.sort(key=lambda x: -x["tpm"])
with open(f"{out}/analysis/expression_top.tsv", "w") as f:
    f.write("transcript\ttpm\treads\tlength\n")
    for q in quant_data[:50]:
        f.write(f"{q['name']}\t{q['tpm']}\t{q['count']}\t{q['length']}\n")

highly_expressed = sum(1 for q in quant_data if q["tpm"] > 1000)
with open(f"{out}/analysis/expression_summary.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"quantified_transcripts\t{len(quant_data)}\n")
    f.write(f"expressed_transcripts\t{expressed}\n")
    f.write(f"highly_expressed\t{highly_expressed}\n")

print(f"  Expression: {expressed}/{len(quant_data)} expressed, {highly_expressed} highly expressed (TPM>1000)")
PYEOF

# 8c: Novel gene detection
python3 << 'PYEOF'
import os

out = os.environ.get("OUT", "outputs")

# Identify novel genes/loci from gffcompare
novel_genes = []
with open(f"{out}/comparison/gffcmp.annotated.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "transcript":
            continue
        attrs = parts[8]
        class_code = ""
        gene_id = ""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("class_code"):
                class_code = attr.split('"')[1] if '"' in attr else ""
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1] if '"' in attr else ""
        if class_code in ("u", "x", "i", "p"):
            novel_genes.append({"gene": gene_id, "class": class_code,
                              "chrom": parts[0], "start": parts[3], "end": parts[4]})

with open(f"{out}/analysis/novel_genes.tsv", "w") as f:
    f.write("gene_id\tclass_code\tchrom\tstart\tend\n")
    for g in novel_genes[:50]:
        f.write(f"{g['gene']}\t{g['class']}\t{g['chrom']}\t{g['start']}\t{g['end']}\n")

print(f"  Novel genes/loci: {len(novel_genes)}")
PYEOF

# ─── Level 9: CONVERGENCE 3 ───
echo "[Level 9 / CONVERGENCE 3] Merging quantification + expression + novel detection..."

# ─── Level 10: CONVERGENCE 4 — Final report ───
echo "[Level 10 / CONVERGENCE 4] Generating final report..."

# Read all analysis files
python3 << PYEOF
import os

out = os.environ.get("OUT", "outputs")
results = os.environ.get("RESULTS", "results")
os.makedirs(results, exist_ok=True)

# Read classification summary
classification = {}
with open(f"{out}/analysis/classification_summary.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        classification[k] = v

# Read transcript stats
tx_stats = {}
with open(f"{out}/analysis/transcript_stats.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        tx_stats[k] = v

# Read expression summary
expr_stats = {}
with open(f"{out}/analysis/expression_summary.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        expr_stats[k] = v

# Count novel genes
novel_count = 0
with open(f"{out}/analysis/novel_genes.tsv") as f:
    next(f)
    for line in f:
        if line.strip():
            novel_count += 1

# Get isoform class distribution
class_dist = {}
with open(f"{out}/analysis/isoform_classification.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            class_dist[parts[1]] = int(parts[2])

# Base-level accuracy from gffcompare
sensitivity = "${SENSITIVITY_BASE}"
precision = "${PRECISION_BASE}"

with open(f"{results}/report.csv", "w") as f:
    f.write("metric,value\n")
    # QC metrics
    f.write(f"total_reads,${TOTAL_READS}\n")
    f.write(f"mean_read_length,${MEAN_LENGTH}\n")
    f.write(f"median_read_length,${MEDIAN_LENGTH}\n")
    f.write(f"read_n50,${N50_READ}\n")
    f.write(f"mean_quality,${MEAN_QUAL}\n")
    f.write(f"total_bases,${TOTAL_BASES}\n")
    # Alignment
    f.write(f"mapped_reads,${MAPPED_READS}\n")
    f.write(f"mapping_pct,${MAPPING_PCT}\n")
    # Assembly
    f.write(f"assembled_transcripts,{tx_stats.get('total_transcripts','0')}\n")
    f.write(f"assembled_genes,{tx_stats.get('total_genes','0')}\n")
    f.write(f"mean_transcript_length,{tx_stats.get('mean_transcript_length','0')}\n")
    f.write(f"transcript_n50,{tx_stats.get('transcript_n50','0')}\n")
    f.write(f"mean_coverage,{tx_stats.get('mean_coverage','0')}\n")
    # Comparison with reference
    f.write(f"base_sensitivity,{sensitivity}\n")
    f.write(f"base_precision,{precision}\n")
    f.write(f"known_isoforms,{classification.get('known_isoforms','0')}\n")
    f.write(f"novel_isoforms,{classification.get('novel_isoforms','0')}\n")
    f.write(f"multi_exon_transcripts,{classification.get('multi_exon_transcripts','0')}\n")
    # Classification
    f.write(f"exact_match_transcripts,{class_dist.get('exact_match',0)}\n")
    f.write(f"novel_junction_transcripts,{class_dist.get('novel_junction',0)}\n")
    f.write(f"intergenic_transcripts,{class_dist.get('intergenic',0)}\n")
    # Expression
    f.write(f"expressed_transcripts,{expr_stats.get('expressed_transcripts','0')}\n")
    f.write(f"highly_expressed_transcripts,{expr_stats.get('highly_expressed','0')}\n")
    # Novel
    f.write(f"novel_loci,{novel_count}\n")

print("Report written to results/report.csv")
PYEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RESULTS}/report.csv"
