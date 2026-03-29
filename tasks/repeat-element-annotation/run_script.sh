#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Task: repeat-element-annotation
# DAG Structure (depth=10, convergence=4, tools=10):
#
#  genome.fasta
#        |
#  +-----+-------------------+
#  |     |                   |
# [RepeatModeler  [RepeatMasker    [python               Level 1-3
#  (de novo TE     (known TE       genome stats
#   library        library)]       (GC, size)]
#   building)]         |
#  |(slow: 1-2h)  [parse RM                              Level 4
#  |               output]
#  |                   |
#  +-------+-----------+
#          |
#  [CONVERGENCE 1]                                       Level 5
#  (de novo + known TE libraries merged)
#  [cat + cd-hit dedup libraries]
#          |
#  [RepeatMasker (comprehensive)]                        Level 6
#  (merged library)
#          |
#  +-------+-------------------+
#  |       |                   |
# [perl   [bedtools          [python                     Level 7
#  parseRM  intersect          TE age
#  landscape] (TE vs genes)]   distribution
#  |       |                   (Kimura)]
#  |       |                   |
#  +-------+-----+-----+------+
#                 |
#         [CONVERGENCE 2]                                Level 8
#         (landscape + gene proximity + age)
#                 |
#         +-------+-----------+
#         |       |           |
#   [python    [python      [python                      Level 9
#    TE class   gene         solo-LTR
#    summary    disruption   detection]
#    (SINE/     analysis]
#     LINE/
#     LTR/DNA)]
#         |       |           |
#         +-------+-----------+
#                 |
#         [CONVERGENCE 3]
#                 |
#         [CONVERGENCE 4] <-- genome stats               Level 10
#         [python report]
# ============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="$WORKDIR/data"
REF="$WORKDIR/reference"
OUT="$WORKDIR/outputs"
RES="$WORKDIR/results"

mkdir -p "$OUT"/{denovo,known,merged,comprehensive,landscape,intersect,age,analysis}
mkdir -p "$RES"

# ============================================================================
# LEVEL 1-3a: De novo TE library construction with RepeatModeler
# ============================================================================
if [ ! -f "$OUT/denovo/genome-families.fa" ]; then
  echo "[Level 1-3a] Running RepeatModeler (de novo TE library)..."
  cd "$OUT/denovo"

  # Build database
  if [ ! -f genome.nsq ] && [ ! -f genome.nhr ]; then
    BuildDatabase -name genome "$DATA/genome.fasta"
  fi

  # Run RepeatModeler
  RepeatModeler -database genome -threads "$THREADS" -LTRStruct 2>&1 | tail -20 || true

  # RepeatModeler outputs to genome-families.fa
  if [ ! -f genome-families.fa ]; then
    # Check for RM output directories
    RMDIR=$(ls -d RM_* 2>/dev/null | tail -1)
    if [ -n "$RMDIR" ] && [ -f "${RMDIR}/consensi.fa.classified" ]; then
      cp "${RMDIR}/consensi.fa.classified" genome-families.fa
    else
      echo "WARNING: RepeatModeler produced no output, creating empty library"
      touch genome-families.fa
    fi
  fi

  cd "$WORKDIR"
fi

# ============================================================================
# LEVEL 1-3b: Known TE library annotation with RepeatMasker
# ============================================================================
if [ ! -f "$OUT/known/genome.fasta.out" ]; then
  echo "[Level 1-3b] Running RepeatMasker (known library)..."
  RepeatMasker -species drosophila \
    -pa "$THREADS" \
    -dir "$OUT/known" \
    -gff -xsmall -a \
    "$DATA/genome.fasta" || true
fi

# ============================================================================
# LEVEL 1-3c: Genome statistics
# ============================================================================
if [ ! -f "$OUT/analysis/genome_stats.tsv" ]; then
  echo "[Level 1-3c] Computing genome statistics..."
  python3 << 'PYEOF'
import os

seqs = {}
current = None
with open("data/genome.fasta") as f:
    for line in f:
        if line.startswith(">"):
            current = line[1:].strip().split()[0]
            seqs[current] = ""
        elif current:
            seqs[current] += line.strip()

total_len = sum(len(s) for s in seqs.values())
gc_count = sum(s.count("G") + s.count("C") + s.count("g") + s.count("c") for s in seqs.values())
n_count = sum(s.count("N") + s.count("n") for s in seqs.values())
gc_frac = gc_count / (total_len - n_count) if (total_len - n_count) > 0 else 0

os.makedirs("outputs/analysis", exist_ok=True)
with open("outputs/analysis/genome_stats.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"total_length\t{total_len}\n")
    f.write(f"num_sequences\t{len(seqs)}\n")
    f.write(f"gc_content\t{gc_frac:.4f}\n")
    f.write(f"n_content\t{n_count / total_len:.4f}\n")
    for name, seq in sorted(seqs.items()):
        f.write(f"seq_{name}_length\t{len(seq)}\n")
print(f"Genome: {total_len:,} bp, {len(seqs)} sequences, GC={gc_frac:.3f}")
PYEOF
fi

# ============================================================================
# LEVEL 4: Parse known RepeatMasker output
# ============================================================================
if [ ! -f "$OUT/known/known_te_summary.tsv" ]; then
  echo "[Level 4] Parsing known RM output..."
  python3 << 'PYEOF'
import os, re

rm_out = "outputs/known/genome.fasta.out"
if not os.path.exists(rm_out):
    print("WARNING: No RepeatMasker output found")
    with open("outputs/known/known_te_summary.tsv", "w") as f:
        f.write("class\tcount\ttotal_bp\n")
else:
    te_classes = {}
    with open(rm_out) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue
            parts = line.split()
            if len(parts) >= 11:
                try:
                    te_class = parts[10] if len(parts) > 10 else "Unknown"
                    begin = int(parts[5])
                    end = int(parts[6])
                    length = end - begin + 1
                    if te_class not in te_classes:
                        te_classes[te_class] = {"count": 0, "total_bp": 0}
                    te_classes[te_class]["count"] += 1
                    te_classes[te_class]["total_bp"] += length
                except (ValueError, IndexError):
                    pass

    with open("outputs/known/known_te_summary.tsv", "w") as f:
        f.write("class\tcount\ttotal_bp\n")
        for cls in sorted(te_classes, key=lambda x: te_classes[x]["total_bp"], reverse=True):
            f.write(f"{cls}\t{te_classes[cls]['count']}\t{te_classes[cls]['total_bp']}\n")

    total_te = sum(v["total_bp"] for v in te_classes.values())
    print(f"Known TEs: {len(te_classes)} classes, {total_te:,} bp total")
PYEOF
fi

# ============================================================================
# LEVEL 5 (CONVERGENCE 1): Merge de novo + known libraries
# ============================================================================
if [ ! -f "$OUT/merged/merged_library.fa" ]; then
  echo "[Level 5] CONVERGENCE 1: Merging TE libraries..."

  # Combine de novo and known libraries
  cat "$OUT/denovo/genome-families.fa" > "$OUT/merged/combined_raw.fa" 2>/dev/null || true

  # Add known library from RepeatMasker's database if available
  RM_LIB=$(find "$OUT/known" -name "*.lib" 2>/dev/null | head -1)
  if [ -n "$RM_LIB" ] && [ -f "$RM_LIB" ]; then
    cat "$RM_LIB" >> "$OUT/merged/combined_raw.fa"
  fi

  # If combined library is empty or very small, just use the de novo library
  if [ ! -s "$OUT/merged/combined_raw.fa" ] || [ "$(wc -c < "$OUT/merged/combined_raw.fa")" -lt 100 ]; then
    echo "WARNING: Combined library is empty/small, using de novo only"
    cp "$OUT/denovo/genome-families.fa" "$OUT/merged/merged_library.fa" 2>/dev/null || touch "$OUT/merged/merged_library.fa"
  else
    # Deduplicate with cd-hit-est
    cd-hit-est -i "$OUT/merged/combined_raw.fa" \
      -o "$OUT/merged/merged_library.fa" \
      -c 0.8 -n 5 -T "$THREADS" -M 4000 || \
      cp "$OUT/merged/combined_raw.fa" "$OUT/merged/merged_library.fa"
  fi

  NSEQS=$(grep -c "^>" "$OUT/merged/merged_library.fa" || true)
  echo "  Merged library: ${NSEQS:-0} sequences"
fi

# ============================================================================
# LEVEL 6: Comprehensive RepeatMasker run with merged library
# ============================================================================
if [ ! -f "$OUT/comprehensive/genome.fasta.out" ]; then
  echo "[Level 6] Running comprehensive RepeatMasker..."

  # Use merged library if it has sequences, otherwise use species library
  NSEQS=$(grep -c "^>" "$OUT/merged/merged_library.fa" 2>/dev/null || true)

  if [ "${NSEQS:-0}" -gt 0 ]; then
    RepeatMasker -lib "$OUT/merged/merged_library.fa" \
      -pa "$THREADS" \
      -dir "$OUT/comprehensive" \
      -gff -xsmall -a \
      "$DATA/genome.fasta" || true
  else
    echo "  Merged library empty, using species library for comprehensive run"
    RepeatMasker -species drosophila \
      -pa "$THREADS" \
      -dir "$OUT/comprehensive" \
      -gff -xsmall -a \
      "$DATA/genome.fasta" || true
  fi
fi

# Use best available RepeatMasker output
if [ -f "$OUT/comprehensive/genome.fasta.out" ]; then
  BEST_RM="$OUT/comprehensive"
elif [ -f "$OUT/known/genome.fasta.out" ]; then
  BEST_RM="$OUT/known"
else
  echo "ERROR: No RepeatMasker output available"
  exit 1
fi
echo "Using RM output from: $BEST_RM"

# ============================================================================
# LEVEL 7a: RepeatMasker landscape (repeat divergence)
# ============================================================================
if [ ! -f "$OUT/landscape/landscape.tsv" ]; then
  echo "[Level 7a] Generating repeat landscape..."
  python3 << 'PYEOF'
import os, re

rm_out = os.environ.get("BEST_RM", "outputs/comprehensive") + "/genome.fasta.out"
if not os.path.exists(rm_out):
    rm_out = "outputs/known/genome.fasta.out"

# Parse divergence values from RM output
divergences = {}  # class -> [div values]
with open(rm_out) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("SW") or line.startswith("score"):
            continue
        parts = line.split()
        if len(parts) >= 11:
            try:
                div = float(parts[1])  # %div from consensus
                te_class = parts[10]
                # Simplify class
                major = te_class.split("/")[0] if "/" in te_class else te_class
                if major not in divergences:
                    divergences[major] = []
                divergences[major].append(div)
            except (ValueError, IndexError):
                pass

os.makedirs("outputs/landscape", exist_ok=True)

# Create binned landscape
bins = list(range(0, 51, 2))  # 0-50% in 2% bins
with open("outputs/landscape/landscape.tsv", "w") as f:
    f.write("bin_start\tbin_end\t" + "\t".join(sorted(divergences.keys())) + "\n")
    for i in range(len(bins)-1):
        row = [str(bins[i]), str(bins[i+1])]
        for cls in sorted(divergences.keys()):
            count = sum(1 for d in divergences[cls] if bins[i] <= d < bins[i+1])
            row.append(str(count))
        f.write("\t".join(row) + "\n")

print(f"Landscape: {len(divergences)} TE classes, {sum(len(v) for v in divergences.values())} elements")
PYEOF
fi

# ============================================================================
# LEVEL 7b: TE-gene intersection
# ============================================================================
if [ ! -f "$OUT/intersect/te_gene_intersect.tsv" ]; then
  echo "[Level 7b] Computing TE-gene intersections..."

  # Convert RM GFF to BED
  RM_GFF=$(ls "$BEST_RM"/*.gff 2>/dev/null | head -1)
  if [ -n "$RM_GFF" ] && [ -f "$RM_GFF" ]; then
    grep -v "^#" "$RM_GFF" | awk -F'\t' 'NF>=9{OFS="\t"; print $1, $4-1, $5, $9, ".", $7}' > "$OUT/intersect/te.bed"

    # Intersect TEs with genes
    bedtools intersect -a "$OUT/intersect/te.bed" -b "$REF/genes.bed" -wa -wb > "$OUT/intersect/te_gene_intersect.tsv" || true
    bedtools closest -a "$OUT/intersect/te.bed" -b "$REF/genes.bed" -d > "$OUT/intersect/te_gene_closest.tsv" || true

    echo "  TE-gene intersections: $(wc -l < "$OUT/intersect/te_gene_intersect.tsv" || echo 0)"
    echo "  TE-gene closest: $(wc -l < "$OUT/intersect/te_gene_closest.tsv" || echo 0)"
  else
    echo "WARNING: No GFF output, creating empty intersection"
    touch "$OUT/intersect/te_gene_intersect.tsv"
    touch "$OUT/intersect/te_gene_closest.tsv"
  fi
fi

# ============================================================================
# LEVEL 7c: TE age distribution (Kimura distance)
# ============================================================================
if [ ! -f "$OUT/age/kimura_distribution.tsv" ]; then
  echo "[Level 7c] Computing TE age distribution (Kimura distances)..."
  python3 << 'PYEOF'
import os, re, math

# Parse .align file for Kimura distances
align_file = None
for d in ["outputs/comprehensive", "outputs/known"]:
    af = os.path.join(d, "genome.fasta.align")
    if os.path.exists(af):
        align_file = af
        break

os.makedirs("outputs/age", exist_ok=True)

if align_file:
    # Parse alignment file for CpG-adjusted Kimura distances
    # RM .align format: blocks separated by blank lines
    kimura_by_class = {}
    current_class = None

    with open(align_file) as f:
        for line in f:
            line = line.strip()
            # Look for lines with divergence info
            if line.startswith("Kimura") or "kimura" in line.lower():
                match = re.search(r'([\d.]+)', line)
                if match and current_class:
                    k = float(match.group(1))
                    major = current_class.split("/")[0]
                    if major not in kimura_by_class:
                        kimura_by_class[major] = []
                    kimura_by_class[major].append(k)
            elif "matching" in line or "Complement" in line.split():
                pass

    # Fallback: use %div from .out file as proxy for Kimura distance
    if not kimura_by_class:
        rm_out = None
        for d in ["outputs/comprehensive", "outputs/known"]:
            of = os.path.join(d, "genome.fasta.out")
            if os.path.exists(of):
                rm_out = of
                break
        if rm_out:
            with open(rm_out) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 11:
                        try:
                            div = float(parts[1])
                            te_class = parts[10]
                            major = te_class.split("/")[0] if "/" in te_class else te_class
                            if major not in kimura_by_class:
                                kimura_by_class[major] = []
                            kimura_by_class[major].append(div)
                        except (ValueError, IndexError):
                            pass

    with open("outputs/age/kimura_distribution.tsv", "w") as f:
        f.write("class\tmedian_divergence\tmean_divergence\tn_elements\n")
        for cls in sorted(kimura_by_class):
            vals = kimura_by_class[cls]
            if vals:
                median = sorted(vals)[len(vals)//2]
                mean = sum(vals) / len(vals)
                f.write(f"{cls}\t{median:.2f}\t{mean:.2f}\t{len(vals)}\n")

    print(f"Kimura analysis: {len(kimura_by_class)} classes")
else:
    print("WARNING: No .align file found")
    with open("outputs/age/kimura_distribution.tsv", "w") as f:
        f.write("class\tmedian_divergence\tmean_divergence\tn_elements\n")
PYEOF
fi

# ============================================================================
# LEVEL 8 (CONVERGENCE 2): landscape + gene proximity + age
# ============================================================================
echo "[Level 8] CONVERGENCE 2: landscape + gene proximity + age ready"

# ============================================================================
# LEVEL 9a: TE class summary (SINE/LINE/LTR/DNA)
# ============================================================================
if [ ! -f "$OUT/analysis/te_class_summary.tsv" ]; then
  echo "[Level 9a] Computing TE class summary..."
  python3 << 'PYEOF'
import os

rm_out = None
for d in ["outputs/comprehensive", "outputs/known"]:
    of = os.path.join(d, "genome.fasta.out")
    if os.path.exists(of):
        rm_out = of
        break

classes = {}  # major_class -> {count, bp}
with open(rm_out) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 11:
            try:
                te_class = parts[10]
                begin = int(parts[5])
                end = int(parts[6])
                length = end - begin + 1

                # Major class
                if "/" in te_class:
                    major = te_class.split("/")[0]
                else:
                    major = te_class

                if major not in classes:
                    classes[major] = {"count": 0, "bp": 0, "subclasses": set()}
                classes[major]["count"] += 1
                classes[major]["bp"] += length
                classes[major]["subclasses"].add(te_class)
            except (ValueError, IndexError):
                pass

os.makedirs("outputs/analysis", exist_ok=True)
with open("outputs/analysis/te_class_summary.tsv", "w") as f:
    f.write("major_class\tcount\ttotal_bp\tnum_subclasses\n")
    for cls in sorted(classes, key=lambda x: classes[x]["bp"], reverse=True):
        f.write(f"{cls}\t{classes[cls]['count']}\t{classes[cls]['bp']}\t{len(classes[cls]['subclasses'])}\n")

total = sum(v["bp"] for v in classes.values())
print(f"TE classes: {len(classes)} major, {total:,} bp total")
PYEOF
fi

# ============================================================================
# LEVEL 9b: Gene disruption analysis
# ============================================================================
if [ ! -f "$OUT/analysis/gene_disruption.tsv" ]; then
  echo "[Level 9b] Analyzing gene disruptions by TEs..."
  python3 << 'PYEOF'
import os

intersect_file = "outputs/intersect/te_gene_intersect.tsv"
closest_file = "outputs/intersect/te_gene_closest.tsv"

gene_disruptions = {}
te_near_genes = {"intronic": 0, "exonic": 0, "upstream_1kb": 0, "upstream_5kb": 0, "intergenic": 0}

if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
    with open(intersect_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 10:
                gene_id = parts[9] if len(parts) > 9 else "unknown"
                if gene_id not in gene_disruptions:
                    gene_disruptions[gene_id] = 0
                gene_disruptions[gene_id] += 1
    te_near_genes["exonic"] = len(gene_disruptions)  # TEs overlapping gene bodies

if os.path.exists(closest_file) and os.path.getsize(closest_file) > 0:
    with open(closest_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 13:
                try:
                    dist = int(parts[-1])
                    if dist == 0:
                        pass  # already counted
                    elif dist <= 1000:
                        te_near_genes["upstream_1kb"] += 1
                    elif dist <= 5000:
                        te_near_genes["upstream_5kb"] += 1
                    else:
                        te_near_genes["intergenic"] += 1
                except ValueError:
                    pass

os.makedirs("outputs/analysis", exist_ok=True)
with open("outputs/analysis/gene_disruption.tsv", "w") as f:
    f.write("category\tcount\n")
    for cat, count in te_near_genes.items():
        f.write(f"{cat}\t{count}\n")
    f.write(f"genes_with_te_insertions\t{len(gene_disruptions)}\n")

print(f"Gene disruptions: {len(gene_disruptions)} genes affected")
PYEOF
fi

# ============================================================================
# LEVEL 9c: Solo-LTR detection
# ============================================================================
if [ ! -f "$OUT/analysis/solo_ltr.tsv" ]; then
  echo "[Level 9c] Detecting solo-LTRs..."
  python3 << 'PYEOF'
import os

rm_out = None
for d in ["outputs/comprehensive", "outputs/known"]:
    of = os.path.join(d, "genome.fasta.out")
    if os.path.exists(of):
        rm_out = of
        break

ltr_elements = []
with open(rm_out) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 11:
            try:
                te_class = parts[10]
                begin = int(parts[5])
                end = int(parts[6])
                length = end - begin + 1
                if "LTR" in te_class:
                    ltr_elements.append({
                        "class": te_class,
                        "length": length,
                        "chrom": parts[4],
                        "start": begin,
                        "end": end
                    })
            except (ValueError, IndexError):
                pass

# Heuristic: solo-LTRs are typically short (<1000 bp) while intact LTR elements are longer
solo_count = sum(1 for e in ltr_elements if e["length"] < 1000)
intact_count = sum(1 for e in ltr_elements if e["length"] >= 1000)

os.makedirs("outputs/analysis", exist_ok=True)
with open("outputs/analysis/solo_ltr.tsv", "w") as f:
    f.write("category\tcount\n")
    f.write(f"total_ltr_elements\t{len(ltr_elements)}\n")
    f.write(f"solo_ltr_candidates\t{solo_count}\n")
    f.write(f"intact_ltr_candidates\t{intact_count}\n")
    if intact_count > 0:
        f.write(f"solo_to_intact_ratio\t{solo_count/intact_count:.2f}\n")

print(f"LTR analysis: {len(ltr_elements)} total, {solo_count} solo, {intact_count} intact")
PYEOF
fi

# ============================================================================
# LEVEL 9-10 (CONVERGENCE 3+4): Final report
# ============================================================================
echo "[Level 10] CONVERGENCE 3+4: Generating final report..."

python3 << 'PYEOF'
import os
import pandas as pd

results = {}

# ---- Genome stats ----
gs = pd.read_csv("outputs/analysis/genome_stats.tsv", sep="\t")
for _, row in gs.iterrows():
    results[row["metric"]] = row["value"]

# ---- TE summary from comprehensive (or known) RM ----
rm_out = None
for d in ["outputs/comprehensive", "outputs/known"]:
    of = os.path.join(d, "genome.fasta.out")
    if os.path.exists(of):
        rm_out = of
        break

total_te_bp = 0
total_te_count = 0
te_classes = set()
if rm_out:
    with open(rm_out) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 11:
                try:
                    begin = int(parts[5])
                    end = int(parts[6])
                    total_te_bp += end - begin + 1
                    total_te_count += 1
                    te_class = parts[10]
                    major = te_class.split("/")[0] if "/" in te_class else te_class
                    te_classes.add(major)
                except (ValueError, IndexError):
                    pass

results["total_repeat_elements"] = total_te_count
results["total_repeat_bp"] = total_te_bp
genome_size = int(results.get("total_length", 1))
results["repeat_fraction"] = round(total_te_bp / genome_size, 4) if genome_size > 0 else 0
results["num_major_te_classes"] = len(te_classes)

# ---- De novo library ----
denovo_lib = "outputs/denovo/genome-families.fa"
if os.path.exists(denovo_lib):
    denovo_count = sum(1 for l in open(denovo_lib) if l.startswith(">"))
    results["denovo_te_families"] = denovo_count

# ---- Merged library ----
merged_lib = "outputs/merged/merged_library.fa"
if os.path.exists(merged_lib):
    merged_count = sum(1 for l in open(merged_lib) if l.startswith(">"))
    results["merged_library_size"] = merged_count

# ---- TE class summary ----
cls_file = "outputs/analysis/te_class_summary.tsv"
if os.path.exists(cls_file):
    cls_df = pd.read_csv(cls_file, sep="\t")
    for _, row in cls_df.iterrows():
        safe_name = str(row["major_class"]).replace("/", "_").replace("-", "_").lower()
        # Only include major TE classes, skip simple repeats
        if safe_name not in ["simple_repeat", "low_complexity", "satellite", "rrna", "trna", "snrna"]:
            results[f"te_{safe_name}_bp"] = int(row["total_bp"])
            results[f"te_{safe_name}_count"] = int(row["count"])

# ---- TE age distribution ----
age_file = "outputs/age/kimura_distribution.tsv"
if os.path.exists(age_file) and os.path.getsize(age_file) > 50:
    age_df = pd.read_csv(age_file, sep="\t")
    if len(age_df) > 0:
        results["te_mean_divergence"] = round(age_df["mean_divergence"].mean(), 2)
        youngest = age_df.loc[age_df["mean_divergence"].idxmin()]
        results["youngest_te_class"] = youngest["class"]

# ---- Gene disruption ----
gene_file = "outputs/analysis/gene_disruption.tsv"
if os.path.exists(gene_file):
    gd = pd.read_csv(gene_file, sep="\t")
    for _, row in gd.iterrows():
        results[f"te_gene_{row['category']}"] = int(row["count"])

# ---- Solo-LTR ----
ltr_file = "outputs/analysis/solo_ltr.tsv"
if os.path.exists(ltr_file):
    ltr = pd.read_csv(ltr_file, sep="\t")
    for _, row in ltr.iterrows():
        results[row["category"]] = row["count"]

# ---- Landscape ----
land_file = "outputs/landscape/landscape.tsv"
if os.path.exists(land_file):
    results["landscape_generated"] = "yes"

# ---- Write CSV ----
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
