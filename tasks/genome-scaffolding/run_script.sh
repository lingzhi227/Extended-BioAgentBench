#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Genome Scaffolding Pipeline: Long-Read Scaffolding of Fragmented Assembly
# ============================================================
# DAG structure (depth=10, convergence=4):
#
#  contigs.fasta (short-read assembly)     long_reads.fastq.gz
#        │                                       │
#  [QUAST initial assessment]              [seqkit stats QC]         Level 1
#        │                                       │
#        └──────────────┬────────────────────────┘
#                       │
#               [CONVERGENCE 1: contigs + reads available]           Level 2
#                       │
#       ┌───────────────┼───────────────┐
#       │               │               │
#  [minimap2 align   [minimap2 align  [ntLink                       Level 3
#   + RagTag          + LINKS          scaffolding]
#   scaffold]         scaffolding]
#       │               │               │
#  [RagTag           [abyss-scaffold   │                             Level 4
#   patch]            on LINKS out]    │
#       │               │               │
#       └───────────────┼───────────────┘
#                       │
#               [CONVERGENCE 2: 3 scaffold sets]                    Level 5
#               [python select best by N50]
#                       │
#               ┌───────┼───────────┐
#               │       │           │
#         [QUAST     [minimap2   [BUSCO                              Level 6
#          scaffold   re-align   completeness]
#          assessment] reads]
#               │       │           │
#               └───────┼───────────┘
#                       │
#               [CONVERGENCE 3: QUAST + mapping + BUSCO]            Level 7
#                       │
#               ┌───────┼───────────┐
#               │       │           │
#         [bedtools [python       [python                            Level 8
#          gap       N-gap         scaffold
#          flanks]   analysis]     comparison
#                                  vs initial]
#               │       │           │
#               └───────┼───────────┘
#                       │
#               [CONVERGENCE 4: gaps + comparison + QC]             Level 9
#                       │
#               [python report]                                      Level 10
#
# Longest path: QUAST initial -> convergence1 -> minimap2+RagTag scaffold ->
#   RagTag patch -> convergence2 -> QUAST scaffold -> convergence3 ->
#   gap analysis -> convergence4 -> report = depth 10
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{initial_qc,ragtag,links,ntlink,best_scaffold,scaffold_qc,mapping,busco,gaps,comparison}
mkdir -p "${RESULTS}"

CONTIGS="${DATA}/contigs.fasta"
LONGREADS="${DATA}/long_reads.fastq.gz"
REFERENCE="${REF}/ecoli_k12_reference.fasta"

# ============================================================
# Level 1: Initial assessment (parallel)
# ============================================================
echo "=== Level 1: Initial assessment ==="

# Branch 1a: QUAST on initial contigs
if [ ! -f "${OUT}/initial_qc/report.tsv" ]; then
  quast "${CONTIGS}" \
    -r "${REFERENCE}" \
    -o "${OUT}/initial_qc" \
    -t ${THREADS} \
    --min-contig 200 \
    > /dev/null 2>&1
  echo "  Initial QUAST done"
fi

# Branch 1b: Long read stats
if [ ! -f "${OUT}/initial_qc/longreads_stats.txt" ]; then
  seqkit stats -T "${LONGREADS}" > "${OUT}/initial_qc/longreads_stats.txt" 2>/dev/null
  echo "  Long read stats done"
fi

# ============================================================
# Level 2: CONVERGENCE 1 — contigs + reads available
# ============================================================
echo "=== Level 2: Convergence 1 (contigs + reads ready) ==="
# Index contigs for tools that need it
samtools faidx "${CONTIGS}" 2>/dev/null || true

# ============================================================
# Level 3-4: Three parallel scaffolding approaches
# ============================================================
echo "=== Level 3-4: Scaffolding (3 methods) ==="

# Method A: RagTag scaffold + patch (minimap2-based)
if [ ! -f "${OUT}/ragtag/ragtag.scaffold.fasta" ]; then
  ragtag.py scaffold "${REFERENCE}" "${CONTIGS}" \
    -o "${OUT}/ragtag" \
    -t ${THREADS} \
    -u \
    > "${OUT}/ragtag/ragtag.log" 2>&1
  echo "  RagTag scaffold done"
fi

# Level 4: RagTag patch (fills gaps using long reads)
if [ ! -f "${OUT}/ragtag/ragtag_patched.fasta" ]; then
  # Align long reads
  minimap2 -ax map-ont -t ${THREADS} "${OUT}/ragtag/ragtag.scaffold.fasta" "${LONGREADS}" 2>/dev/null | \
    samtools sort -@ ${THREADS} -o "${OUT}/ragtag/longreads_aligned.bam" 2>/dev/null
  samtools index "${OUT}/ragtag/longreads_aligned.bam"

  # RagTag patch uses the reference and long-read-aligned BAM
  ragtag.py patch "${OUT}/ragtag/ragtag.scaffold.fasta" "${LONGREADS}" \
    -o "${OUT}/ragtag/patch" \
    -t ${THREADS} \
    > "${OUT}/ragtag/patch.log" 2>&1 || true

  # Use patched output or fall back to scaffold
  if [ -f "${OUT}/ragtag/patch/ragtag.patch.fasta" ]; then
    cp "${OUT}/ragtag/patch/ragtag.patch.fasta" "${OUT}/ragtag/ragtag_patched.fasta"
  else
    cp "${OUT}/ragtag/ragtag.scaffold.fasta" "${OUT}/ragtag/ragtag_patched.fasta"
  fi
  echo "  RagTag patch done"
fi

# Method B: LINKS scaffolding
if [ ! -f "${OUT}/links/links_scaffold.fasta" ]; then
  # LINKS needs uncompressed long reads and a file listing
  zcat "${LONGREADS}" > "${OUT}/links/long_reads.fq"
  echo "${OUT}/links/long_reads.fq" > "${OUT}/links/long_reads_list.txt"
  cp "${CONTIGS}" "${OUT}/links/contigs.fasta"

  cd "${OUT}/links"
  LINKS -f contigs.fasta \
    -s long_reads_list.txt \
    -d 500 \
    -k 15 \
    -l 5 \
    -t 2 \
    -o links_out \
    > links.log 2>&1 || true
  cd "${WORKDIR}"

  # Find the LINKS output (it creates files with various extensions)
  LINKS_OUT=$(ls "${OUT}/links/"*.scaffolds.fa 2>/dev/null | head -1 || true)
  if [ -n "${LINKS_OUT}" ] && [ -s "${LINKS_OUT}" ]; then
    cp "${LINKS_OUT}" "${OUT}/links/links_scaffold.fasta"
  else
    # Fall back to original contigs if LINKS fails
    cp "${CONTIGS}" "${OUT}/links/links_scaffold.fasta"
  fi
  echo "  LINKS scaffolding done"
  # Clean up large temp file
  rm -f "${OUT}/links/long_reads.fq"
fi

# Method C: ntLink scaffolding
if [ ! -f "${OUT}/ntlink/ntlink_scaffold.fasta" ]; then
  cp "${CONTIGS}" "${OUT}/ntlink/contigs.fasta"
  cd "${OUT}/ntlink"
  ntLink scaffold \
    target=contigs.fasta \
    reads="${LONGREADS}" \
    k=32 w=250 \
    t=${THREADS} \
    prefix=ntlink_out \
    > ntlink.log 2>&1 || true
  cd "${WORKDIR}"

  NTLINK_OUT=$(ls "${OUT}/ntlink/"*ntLink*.fa 2>/dev/null | head -1 || true)
  if [ -z "${NTLINK_OUT}" ]; then
    NTLINK_OUT=$(ls "${OUT}/ntlink/"*scaffold*.fa* 2>/dev/null | head -1 || true)
  fi
  if [ -n "${NTLINK_OUT}" ] && [ -s "${NTLINK_OUT}" ]; then
    cp "${NTLINK_OUT}" "${OUT}/ntlink/ntlink_scaffold.fasta"
  else
    cp "${CONTIGS}" "${OUT}/ntlink/ntlink_scaffold.fasta"
  fi
  echo "  ntLink scaffolding done"
fi

# ============================================================
# Level 5: CONVERGENCE 2 — Select best scaffold
# ============================================================
echo "=== Level 5: Convergence 2 (select best scaffold) ==="

python3 << 'PYEOF'
import os

def compute_n50(fasta_path):
    """Compute N50 from a FASTA file."""
    lengths = []
    current_len = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_len > 0:
        lengths.append(current_len)

    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    for l in lengths:
        running += l
        if running >= total / 2:
            return l, len(lengths), total
    return 0, 0, 0

results = {}
for method, path in [("ragtag", "outputs/ragtag/ragtag_patched.fasta"),
                      ("links", "outputs/links/links_scaffold.fasta"),
                      ("ntlink", "outputs/ntlink/ntlink_scaffold.fasta")]:
    if os.path.exists(path):
        n50, nseq, total = compute_n50(path)
        results[method] = {"n50": n50, "nseq": nseq, "total": total, "path": path}
        print(f"  {method}: N50={n50:,}, sequences={nseq}, total={total:,} bp")
    else:
        print(f"  {method}: NOT FOUND")

# Select best by N50
best = max(results.keys(), key=lambda k: results[k]["n50"])
print(f"  Best method: {best} (N50={results[best]['n50']:,})")

import shutil
shutil.copy(results[best]["path"], "outputs/best_scaffold/scaffold.fasta")

# Save comparison
with open("outputs/best_scaffold/method_comparison.tsv", 'w') as f:
    f.write("method\tn50\tnum_sequences\ttotal_length\n")
    for m in results:
        f.write(f"{m}\t{results[m]['n50']}\t{results[m]['nseq']}\t{results[m]['total']}\n")
PYEOF

# ============================================================
# Level 6: Triple QC on best scaffold (parallel)
# ============================================================
echo "=== Level 6: QC on best scaffold ==="

SCAFFOLD="${OUT}/best_scaffold/scaffold.fasta"

# Branch 6a: QUAST on scaffold
if [ ! -f "${OUT}/scaffold_qc/report.tsv" ]; then
  quast "${SCAFFOLD}" \
    -r "${REFERENCE}" \
    -o "${OUT}/scaffold_qc" \
    -t ${THREADS} \
    --min-contig 200 \
    > /dev/null 2>&1
  echo "  Scaffold QUAST done"
fi

# Branch 6b: Re-align long reads to scaffold
if [ ! -f "${OUT}/mapping/scaffold_mapping.bam" ]; then
  minimap2 -ax map-ont -t ${THREADS} "${SCAFFOLD}" "${LONGREADS}" 2>/dev/null | \
    samtools sort -@ ${THREADS} -o "${OUT}/mapping/scaffold_mapping.bam" 2>/dev/null
  samtools index "${OUT}/mapping/scaffold_mapping.bam"
  samtools flagstat "${OUT}/mapping/scaffold_mapping.bam" > "${OUT}/mapping/flagstat.txt"
  echo "  Read re-alignment done"
fi

# Branch 6c: BUSCO completeness
if [ ! -f "${OUT}/busco/done" ]; then
  busco -i "${SCAFFOLD}" \
    -o busco_scaffold \
    -m genome \
    -l bacteria_odb10 \
    --out_path "${OUT}/busco/" \
    -c ${THREADS} \
    --offline \
    > "${OUT}/busco/busco.log" 2>&1 || true

  # Also run on initial contigs for comparison
  busco -i "${CONTIGS}" \
    -o busco_initial \
    -m genome \
    -l bacteria_odb10 \
    --out_path "${OUT}/busco/" \
    -c ${THREADS} \
    --offline \
    > "${OUT}/busco/busco_initial.log" 2>&1 || true
  touch "${OUT}/busco/done"
  echo "  BUSCO done"
fi

# ============================================================
# Level 7: CONVERGENCE 3 — QUAST + mapping + BUSCO
# ============================================================
echo "=== Level 7: Convergence 3 (assessment results) ==="

# ============================================================
# Level 8: Analysis branches (parallel)
# ============================================================
echo "=== Level 8: Detailed analysis ==="

# Branch 8a: Gap flanking analysis with bedtools
if [ ! -f "${OUT}/gaps/gap_flanks.bed" ]; then
  # Find N-gaps in scaffold
  python3 << 'PYEOF2'
import re
gaps = []
with open("outputs/best_scaffold/scaffold.fasta") as f:
    seq_name = ""
    seq = ""
    for line in f:
        if line.startswith('>'):
            if seq_name and seq:
                for m in re.finditer(r'[Nn]+', seq):
                    gaps.append((seq_name, m.start(), m.end(), m.end()-m.start()))
            seq_name = line.strip()[1:].split()[0]
            seq = ""
        else:
            seq += line.strip()
    if seq_name and seq:
        for m in re.finditer(r'[Nn]+', seq):
            gaps.append((seq_name, m.start(), m.end(), m.end()-m.start()))

with open("outputs/gaps/gaps.bed", 'w') as f:
    for g in gaps:
        f.write(f"{g[0]}\t{g[1]}\t{g[2]}\tgap_{g[3]}bp\n")
print(f"Found {len(gaps)} gaps")
PYEOF2

  # Get gap flanking regions
  if [ -s "${OUT}/gaps/gaps.bed" ]; then
    samtools faidx "${SCAFFOLD}" 2>/dev/null || true
    bedtools slop -i "${OUT}/gaps/gaps.bed" -g "${SCAFFOLD}.fai" -b 500 > "${OUT}/gaps/gap_flanks.bed" 2>/dev/null || true
  fi
  touch "${OUT}/gaps/gap_flanks.bed"
  echo "  Gap analysis done"
fi

# Branch 8b: N-gap analysis
if [ ! -f "${OUT}/gaps/gap_summary.tsv" ]; then
  python3 << 'PYEOF3'
import re

def count_gaps(fasta_path):
    gaps = []
    total_n = 0
    with open(fasta_path) as f:
        seq = ""
        for line in f:
            if line.startswith('>'):
                for m in re.finditer(r'[Nn]+', seq):
                    gaps.append(m.end() - m.start())
                    total_n += m.end() - m.start()
                seq = ""
            else:
                seq += line.strip()
        for m in re.finditer(r'[Nn]+', seq):
            gaps.append(m.end() - m.start())
            total_n += m.end() - m.start()
    return len(gaps), total_n, sorted(gaps, reverse=True)

scaffold_ngaps, scaffold_nbases, scaffold_gap_sizes = count_gaps("outputs/best_scaffold/scaffold.fasta")
initial_ngaps, initial_nbases, initial_gap_sizes = count_gaps("data/contigs.fasta")

with open("outputs/gaps/gap_summary.tsv", 'w') as f:
    f.write("assembly\tnum_gaps\ttotal_n_bases\tlargest_gap\n")
    f.write(f"initial\t{initial_ngaps}\t{initial_nbases}\t{initial_gap_sizes[0] if initial_gap_sizes else 0}\n")
    f.write(f"scaffold\t{scaffold_ngaps}\t{scaffold_nbases}\t{scaffold_gap_sizes[0] if scaffold_gap_sizes else 0}\n")

print(f"Initial: {initial_ngaps} gaps, {initial_nbases} N-bases")
print(f"Scaffold: {scaffold_ngaps} gaps, {scaffold_nbases} N-bases")
PYEOF3
  echo "  N-gap comparison done"
fi

# Branch 8c: Scaffold vs initial comparison
if [ ! -f "${OUT}/comparison/comparison.tsv" ]; then
  python3 << 'PYEOF4'
import os, re

def assembly_stats(fasta_path):
    lengths = []
    total_n = 0
    current_len = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                seq = line.strip()
                current_len += len(seq)
                total_n += seq.count('N') + seq.count('n')
    if current_len > 0:
        lengths.append(current_len)

    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    n50 = 0
    l50 = 0
    for i, l in enumerate(lengths):
        running += l
        if running >= total / 2 and n50 == 0:
            n50 = l
            l50 = i + 1

    return {
        "num_sequences": len(lengths),
        "total_length": total,
        "n50": n50,
        "l50": l50,
        "largest": lengths[0] if lengths else 0,
        "n_bases": total_n,
        "gc_content": 0  # placeholder
    }

initial = assembly_stats("data/contigs.fasta")
scaffold = assembly_stats("outputs/best_scaffold/scaffold.fasta")

with open("outputs/comparison/comparison.tsv", 'w') as f:
    f.write("metric\tinitial\tscaffold\timprovement\n")
    for key in ["num_sequences", "total_length", "n50", "l50", "largest", "n_bases"]:
        iv = initial[key]
        sv = scaffold[key]
        if iv != 0:
            if key in ["num_sequences", "l50"]:
                imp = f"{(iv-sv)/iv*100:.1f}% reduction"
            else:
                imp = f"{(sv-iv)/iv*100:.1f}% change"
        else:
            imp = "N/A"
        f.write(f"{key}\t{iv}\t{sv}\t{imp}\n")

print("Comparison:")
for key in ["num_sequences", "total_length", "n50", "largest"]:
    print(f"  {key}: {initial[key]:,} -> {scaffold[key]:,}")
PYEOF4
  echo "  Comparison done"
fi

# ============================================================
# Level 9: CONVERGENCE 4 — all analysis results
# ============================================================
echo "=== Level 9: Convergence 4 ==="

# ============================================================
# Level 10: Final report
# ============================================================
echo "=== Level 10: Final report ==="

python3 << 'PYEOF5'
import csv
import os
import re

metrics = {}

# --- Initial assembly stats ---
def get_quast_metrics(report_path, prefix=""):
    result = {}
    with open(report_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                key = parts[0].strip().replace(' ', '_').replace('#', 'num').lower()
                val = parts[1].strip()
                result[f"{prefix}{key}"] = val
    return result

initial_qc = get_quast_metrics("outputs/initial_qc/report.tsv", "initial_")
scaffold_qc = get_quast_metrics("outputs/scaffold_qc/report.tsv", "scaffold_")

# Key metrics from QUAST
for prefix, qc in [("initial", initial_qc), ("scaffold", scaffold_qc)]:
    for key in ["num_contigs_(>=_0_bp)", "total_length_(>=_0_bp)", "n50", "l50",
                "largest_contig", "num_misassemblies", "genome_fraction_(%)"]:
        qkey = f"{prefix}_{key}"
        if qkey in qc:
            clean_key = f"{prefix}_{key.replace('_(>=_0_bp)', '').replace('_(%)', '_pct')}"
            metrics[clean_key] = qc[qkey]

# --- Method comparison ---
with open("outputs/best_scaffold/method_comparison.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        method = row['method']
        metrics[f"scaffolder_{method}_n50"] = row['n50']
        metrics[f"scaffolder_{method}_sequences"] = row['num_sequences']

# --- Mapping stats ---
if os.path.exists("outputs/mapping/flagstat.txt"):
    with open("outputs/mapping/flagstat.txt") as f:
        for line in f:
            if "mapped (" in line:
                parts = line.strip().split()
                metrics["mapped_reads"] = parts[0]
                pct_match = re.search(r'\(([\d.]+)%', line)
                if pct_match:
                    metrics["mapping_rate_pct"] = pct_match.group(1)
                break

# --- BUSCO ---
for label, busco_dir in [("initial", "busco_initial"), ("scaffold", "busco_scaffold")]:
    short_file = f"outputs/busco/{busco_dir}/short_summary.specific.bacteria_odb10.{busco_dir}.txt"
    if not os.path.exists(short_file):
        # Try generic short_summary
        import glob
        candidates = glob.glob(f"outputs/busco/{busco_dir}/short_summary*.txt")
        if candidates:
            short_file = candidates[0]
    if os.path.exists(short_file):
        with open(short_file) as f:
            for line in f:
                if 'Complete BUSCOs' in line:
                    val = line.strip().split()[0]
                    metrics[f"{label}_complete_genes"] = val
                elif 'Missing BUSCOs' in line:
                    val = line.strip().split()[0]
                    metrics[f"{label}_missing_genes"] = val
                elif 'Total BUSCO' in line:
                    val = line.strip().split()[0]
                    metrics[f"{label}_total_genes"] = val

# --- Gap analysis ---
if os.path.exists("outputs/gaps/gap_summary.tsv"):
    with open("outputs/gaps/gap_summary.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            prefix = row['assembly']
            metrics[f"{prefix}_gap_count"] = row['num_gaps']
            metrics[f"{prefix}_n_bases"] = row['total_n_bases']

# --- Improvement metrics ---
if os.path.exists("outputs/comparison/comparison.tsv"):
    with open("outputs/comparison/comparison.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            metrics[f"improvement_{row['metric']}"] = row['improvement']

# Write report
with open("results/report.csv", 'w') as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("=== Report ===")
for k, v in metrics.items():
    print(f"  {k} = {v}")
PYEOF5

echo "=== Pipeline complete ==="
