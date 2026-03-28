#!/usr/bin/env bash
set -euo pipefail

# CRISPR Screen Analysis Pipeline
# Data: Human Brunello library screen (APR-246 drug sensitivity)
# Samples: T0 (baseline), T8_Drug (APR-246 treated), T8_Vehicle (DMSO control)
#
# DAG Structure (depth=10, convergence=4):
#
#   [T0.fq]──────────[Drug.fq]──────────[Vehicle.fq]
#     │                  │                    │
#   fastqc             fastqc              fastqc          (Step 1: QC)
#     │                  │                    │
#   cutadapt           cutadapt            cutadapt        (Step 2: Trim)
#     │                  │                    │
#     └──────────────────┼────────────────────┘
#                        │
#                   mageck count                           (Step 3: CONVERGE #1)
#                  ╱       │       ╲
#          mageck test  mageck test  mageck mle            (Step 4: 3-way parallel)
#          (drug/T0)    (veh/T0)    (all conditions)
#                  ╲       │       ╱
#                   merge rankings                         (Step 5: CONVERGE #2)
#                     ╱         ╲
#             drug-specific    count QC metrics             (Step 6: parallel)
#             hit analysis     + gini index
#                     ╲         ╱
#                   pathway enrichment                     (Step 7: CONVERGE #3)
#                        │
#                   hit classification                     (Step 8)
#                        │
#                   final report                           (Step 9: CONVERGE #4 w/ multiqc)
#                        │
#                    report.csv                            (Step 10)

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}/fastqc_raw" "${OUT}/fastqc_trimmed" "${OUT}/trimmed"
mkdir -p "${OUT}/count" "${OUT}/rra_drug" "${OUT}/rra_vehicle" "${OUT}/mle"
mkdir -p "${OUT}/comparison" "${OUT}/multiqc" "${RESULTS}"

LIBRARY="${REF}/library.tsv"

# Vector backbone containing the sgRNA cassette
# Reads have: [variable prefix]ACCG[20bp sgRNA]GTTT[scaffold]
# 5' trim sequence for cutadapt (linked adapter)
VECTOR_5="CTTGTGGAAAGGACGAAACACCG"
SCAFFOLD_3="GTTTTAGAGCTAGAAATAGCAAGTT"

# ============================================================================
# Step 1: FastQC on raw reads (parallel per sample)
# ============================================================================
echo "[Step 1] Running FastQC on raw reads..."
for fq in "${DATA}"/T*.fastq.gz; do
    BASENAME=$(basename "$fq" .fastq.gz)
    if [ ! -f "${OUT}/fastqc_raw/${BASENAME}_fastqc.html" ]; then
        fastqc -t "${THREADS}" -o "${OUT}/fastqc_raw" "$fq" &
    fi
done
wait
echo "[Step 1] FastQC raw done."

# ============================================================================
# Step 2: Cutadapt — extract sgRNA from vector context (parallel per sample)
# Uses linked adapter: trims 5' vector then 3' scaffold, keeping just the sgRNA
# ============================================================================
echo "[Step 2] Extracting sgRNA sequences with cutadapt..."
for fq in "${DATA}"/T*.fastq.gz; do
    BASENAME=$(basename "$fq" .fastq.gz)
    TRIMMED="${OUT}/trimmed/${BASENAME}_trimmed.fastq.gz"
    if [ ! -f "$TRIMMED" ]; then
        cutadapt \
            -g "${VECTOR_5}...${SCAFFOLD_3}" \
            -e 0.15 \
            --discard-untrimmed \
            --minimum-length 18 \
            --maximum-length 24 \
            -o "$TRIMMED" \
            "$fq" \
            > "${OUT}/trimmed/${BASENAME}_cutadapt.log" 2>&1 &
    fi
done
wait
echo "[Step 2] Trimming done."

# ============================================================================
# Step 2b: FastQC on trimmed reads
# ============================================================================
echo "[Step 2b] Running FastQC on trimmed reads..."
for fq in "${OUT}/trimmed"/*_trimmed.fastq.gz; do
    BASENAME=$(basename "$fq" .fastq.gz)
    if [ ! -f "${OUT}/fastqc_trimmed/${BASENAME}_fastqc.html" ]; then
        fastqc -t "${THREADS}" -o "${OUT}/fastqc_trimmed" "$fq" &
    fi
done
wait
echo "[Step 2b] FastQC trimmed done."

# ============================================================================
# Step 3: MAGeCK count — CONVERGE all 3 samples (CONVERGENCE #1)
# ============================================================================
echo "[Step 3] Running MAGeCK count (convergence #1: all samples)..."
if [ ! -f "${OUT}/count/screen.count.txt" ]; then
    mageck count \
        -l "$LIBRARY" \
        -n "${OUT}/count/screen" \
        --sample-label "T0,Drug,Vehicle" \
        --fastq \
            "${OUT}/trimmed/T0_control_trimmed.fastq.gz" \
            "${OUT}/trimmed/T8_drug_trimmed.fastq.gz" \
            "${OUT}/trimmed/T8_vehicle_trimmed.fastq.gz" \
        --norm-method median \
        2>&1 | tee "${OUT}/count/mageck_count.log"
fi
echo "[Step 3] MAGeCK count done."

# ============================================================================
# Step 4a: MAGeCK test (RRA) — Drug vs T0
# ============================================================================
echo "[Step 4a] Running MAGeCK test (RRA): Drug vs T0..."
if [ ! -f "${OUT}/rra_drug/drug_vs_t0.gene_summary.txt" ]; then
    mageck test \
        -k "${OUT}/count/screen.count.txt" \
        -t Drug \
        -c T0 \
        -n "${OUT}/rra_drug/drug_vs_t0" \
        --gene-lfc-method alphamedian \
        2>&1 | tee "${OUT}/rra_drug/mageck_test_drug.log"
fi

# ============================================================================
# Step 4b: MAGeCK test (RRA) — Vehicle vs T0
# ============================================================================
echo "[Step 4b] Running MAGeCK test (RRA): Vehicle vs T0..."
if [ ! -f "${OUT}/rra_vehicle/vehicle_vs_t0.gene_summary.txt" ]; then
    mageck test \
        -k "${OUT}/count/screen.count.txt" \
        -t Vehicle \
        -c T0 \
        -n "${OUT}/rra_vehicle/vehicle_vs_t0" \
        --gene-lfc-method alphamedian \
        2>&1 | tee "${OUT}/rra_vehicle/mageck_test_vehicle.log"
fi

# ============================================================================
# Step 4c: MAGeCK MLE — multi-condition modeling
# ============================================================================
echo "[Step 4c] Running MAGeCK MLE..."
# Create design matrix for MLE
DESIGN="${OUT}/mle/design_matrix.txt"
if [ ! -f "$DESIGN" ]; then
    cat > "$DESIGN" << 'DESIGN_EOF'
Samples	baseline	drug	vehicle
T0	1	0	0
Drug	1	1	0
Vehicle	1	0	1
DESIGN_EOF
fi

if [ ! -f "${OUT}/mle/screen_mle.gene_summary.txt" ]; then
    mageck mle \
        -k "${OUT}/count/screen.count.txt" \
        -d "$DESIGN" \
        -n "${OUT}/mle/screen_mle" \
        2>&1 | tee "${OUT}/mle/mageck_mle.log"
fi
echo "[Step 4] All three analysis methods done."

# ============================================================================
# Step 5: Merge rankings — CONVERGENCE #2 (multi-method)
# ============================================================================
echo "[Step 5] Merging rankings from RRA and MLE (convergence #2)..."
python3 << 'MERGE_PY'
import csv
import os

OUT = os.environ.get("OUT", "outputs")

# Load RRA drug results
rra_drug_genes = {}
with open(f"{OUT}/rra_drug/drug_vs_t0.gene_summary.txt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        rra_drug_genes[row['id']] = {
            'neg_rank': int(row['neg|rank']),
            'neg_fdr': float(row['neg|fdr']),
            'neg_lfc': float(row['neg|lfc']),
            'pos_rank': int(row['pos|rank']),
            'pos_fdr': float(row['pos|fdr']),
            'pos_lfc': float(row['pos|lfc']),
        }

# Load RRA vehicle results
rra_veh_genes = {}
with open(f"{OUT}/rra_vehicle/vehicle_vs_t0.gene_summary.txt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        rra_veh_genes[row['id']] = {
            'neg_rank': int(row['neg|rank']),
            'neg_fdr': float(row['neg|fdr']),
            'neg_lfc': float(row['neg|lfc']),
        }

# Load MLE results
mle_genes = {}
mle_file = f"{OUT}/mle/screen_mle.gene_summary.txt"
if os.path.exists(mle_file):
    with open(mle_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row['Gene']
            # MLE has drug|beta, drug|fdr columns
            try:
                mle_genes[gene] = {
                    'drug_beta': float(row.get('drug|beta', 0)),
                    'drug_fdr': float(row.get('drug|fdr', 1)),
                }
            except (ValueError, KeyError):
                pass

# Merge: identify concordant hits
concordance = {}
for gene in rra_drug_genes:
    rra_neg_fdr = rra_drug_genes[gene]['neg_fdr']
    rra_neg_rank = rra_drug_genes[gene]['neg_rank']
    mle_fdr = mle_genes.get(gene, {}).get('drug_fdr', 1.0)
    mle_beta = mle_genes.get(gene, {}).get('drug_beta', 0.0)
    veh_neg_fdr = rra_veh_genes.get(gene, {}).get('neg_fdr', 1.0)
    concordance[gene] = {
        'rra_neg_rank': rra_neg_rank,
        'rra_neg_fdr': rra_neg_fdr,
        'rra_neg_lfc': rra_drug_genes[gene]['neg_lfc'],
        'mle_beta': mle_beta,
        'mle_fdr': mle_fdr,
        'veh_neg_fdr': veh_neg_fdr,
        'drug_specific': rra_neg_fdr < 0.25 and veh_neg_fdr >= 0.25,
        'concordant': rra_neg_fdr < 0.25 and mle_fdr < 0.25,
    }

# Write merged results
os.makedirs(f"{OUT}/comparison", exist_ok=True)
with open(f"{OUT}/comparison/merged_rankings.tsv", 'w') as f:
    f.write("gene\trra_neg_rank\trra_neg_fdr\trra_neg_lfc\tmle_beta\tmle_fdr\tveh_neg_fdr\tdrug_specific\tconcordant\n")
    for gene, data in sorted(concordance.items(), key=lambda x: x[1]['rra_neg_rank']):
        f.write(f"{gene}\t{data['rra_neg_rank']}\t{data['rra_neg_fdr']:.6f}\t{data['rra_neg_lfc']:.4f}\t"
                f"{data['mle_beta']:.4f}\t{data['mle_fdr']:.6f}\t{data['veh_neg_fdr']:.6f}\t"
                f"{data['drug_specific']}\t{data['concordant']}\n")

print(f"Merged {len(concordance)} genes")
drug_specific = sum(1 for g in concordance.values() if g['drug_specific'])
concordant = sum(1 for g in concordance.values() if g['concordant'])
print(f"Drug-specific: {drug_specific}, Concordant (RRA+MLE): {concordant}")
MERGE_PY
echo "[Step 5] Ranking merge done."

# ============================================================================
# Step 6a: Drug-specific hit analysis (parallel)
# ============================================================================
echo "[Step 6a] Identifying drug-specific hits..."
python3 << 'DRUG_SPEC_PY'
import csv
import os

OUT = os.environ.get("OUT", "outputs")

# Read merged rankings
drug_specific_genes = []
with open(f"{OUT}/comparison/merged_rankings.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['drug_specific'] == 'True':
            drug_specific_genes.append({
                'gene': row['gene'],
                'rra_neg_fdr': float(row['rra_neg_fdr']),
                'rra_neg_lfc': float(row['rra_neg_lfc']),
                'veh_neg_fdr': float(row['veh_neg_fdr']),
            })

drug_specific_genes.sort(key=lambda x: x['rra_neg_fdr'])

with open(f"{OUT}/comparison/drug_specific_hits.tsv", 'w') as f:
    f.write("gene\trra_neg_fdr\trra_neg_lfc\tveh_neg_fdr\n")
    for g in drug_specific_genes:
        f.write(f"{g['gene']}\t{g['rra_neg_fdr']:.6f}\t{g['rra_neg_lfc']:.4f}\t{g['veh_neg_fdr']:.6f}\n")

print(f"Drug-specific hits: {len(drug_specific_genes)}")
DRUG_SPEC_PY

# ============================================================================
# Step 6b: Count QC metrics and Gini index (parallel)
# ============================================================================
echo "[Step 6b] Computing count QC metrics..."
python3 << 'QC_PY'
import csv
import math
import os

OUT = os.environ.get("OUT", "outputs")

# Read count table
counts = {"T0": [], "Drug": [], "Vehicle": []}
total_sgrnas = 0
with open(f"{OUT}/count/screen.count.txt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        total_sgrnas += 1
        counts["T0"].append(int(row["T0"]))
        counts["Drug"].append(int(row["Drug"]))
        counts["Vehicle"].append(int(row["Vehicle"]))

def gini(values):
    """Compute Gini index of a distribution."""
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    if n == 0 or sum(sorted_vals) == 0:
        return 0.0
    cumsum = 0
    total = sum(sorted_vals)
    gini_sum = 0
    for i, v in enumerate(sorted_vals):
        cumsum += v
        gini_sum += (2 * (i + 1) - n - 1) * v
    return gini_sum / (n * total)

qc = {}
for sample in ["T0", "Drug", "Vehicle"]:
    vals = counts[sample]
    qc[f"total_counts_{sample.lower()}"] = sum(vals)
    qc[f"zero_count_sgrnas_{sample.lower()}"] = sum(1 for v in vals if v == 0)
    qc[f"gini_index_{sample.lower()}"] = round(gini(vals), 4)
    qc[f"median_count_{sample.lower()}"] = sorted(vals)[len(vals) // 2]

qc["total_sgrnas"] = total_sgrnas

with open(f"{OUT}/comparison/count_qc.tsv", 'w') as f:
    for k, v in qc.items():
        f.write(f"{k}\t{v}\n")
    print(f"  {k}: {v}")

print(f"QC metrics computed for {total_sgrnas} sgRNAs")
QC_PY

# ============================================================================
# Step 7: Pathway enrichment — CONVERGENCE #3 (drug-specific + QC)
# ============================================================================
echo "[Step 7] Running pathway enrichment (convergence #3)..."
python3 << 'PATHWAY_PY'
import csv
import os
from collections import defaultdict

OUT = os.environ.get("OUT", "outputs")

# Load gene rankings
gene_ranks = {}
with open(f"{OUT}/comparison/merged_rankings.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene_ranks[row['gene']] = {
            'rank': int(row['rra_neg_rank']),
            'fdr': float(row['rra_neg_fdr']),
            'lfc': float(row['rra_neg_lfc']),
        }

# Simple enrichment: categorize genes by essentiality based on screen results
categories = {
    'essential_drug': [],  # Depleted in drug (FDR < 0.25)
    'essential_common': [],  # Depleted in both drug and vehicle
    'enriched_drug': [],  # Enriched in drug
    'neutral': [],
}

with open(f"{OUT}/comparison/merged_rankings.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row['gene']
        drug_fdr = float(row['rra_neg_fdr'])
        veh_fdr = float(row['veh_neg_fdr'])

        if drug_fdr < 0.25 and veh_fdr < 0.25:
            categories['essential_common'].append(gene)
        elif drug_fdr < 0.25:
            categories['essential_drug'].append(gene)
        else:
            categories['neutral'].append(gene)

# Write pathway/category summary
with open(f"{OUT}/comparison/gene_categories.tsv", 'w') as f:
    f.write("category\tcount\ttop_genes\n")
    for cat, genes in categories.items():
        top = ','.join(genes[:5]) if genes else 'none'
        f.write(f"{cat}\t{len(genes)}\t{top}\n")
        print(f"  {cat}: {len(genes)} genes")

print("Gene categorization done.")
PATHWAY_PY

# ============================================================================
# Step 8: Hit classification with multi-method consensus
# ============================================================================
echo "[Step 8] Classifying hits with multi-method consensus..."
python3 << 'CLASSIFY_PY'
import csv
import os

OUT = os.environ.get("OUT", "outputs")

# Read merged rankings
hits = []
with open(f"{OUT}/comparison/merged_rankings.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row['gene']
        rra_fdr = float(row['rra_neg_fdr'])
        mle_fdr = float(row['mle_fdr'])
        rra_lfc = float(row['rra_neg_lfc'])
        veh_fdr = float(row['veh_neg_fdr'])

        # Classification tiers
        if rra_fdr < 0.05 and mle_fdr < 0.05:
            tier = "high_confidence"
        elif rra_fdr < 0.25 or mle_fdr < 0.25:
            tier = "moderate_confidence"
        else:
            tier = "not_significant"

        hits.append({
            'gene': gene,
            'rra_fdr': rra_fdr,
            'mle_fdr': mle_fdr,
            'rra_lfc': rra_lfc,
            'tier': tier,
            'drug_specific': rra_fdr < 0.25 and veh_fdr >= 0.25,
        })

# Sort by tier then RRA FDR
tier_order = {'high_confidence': 0, 'moderate_confidence': 1, 'not_significant': 2}
hits.sort(key=lambda x: (tier_order[x['tier']], x['rra_fdr']))

with open(f"{OUT}/comparison/classified_hits.tsv", 'w') as f:
    f.write("gene\ttier\trra_neg_fdr\tmle_fdr\trra_neg_lfc\tdrug_specific\n")
    for h in hits:
        f.write(f"{h['gene']}\t{h['tier']}\t{h['rra_fdr']:.6f}\t{h['mle_fdr']:.6f}\t"
                f"{h['rra_lfc']:.4f}\t{h['drug_specific']}\n")

high = sum(1 for h in hits if h['tier'] == 'high_confidence')
moderate = sum(1 for h in hits if h['tier'] == 'moderate_confidence')
drug_spec = sum(1 for h in hits if h['drug_specific'])
print(f"High confidence: {high}, Moderate: {moderate}, Drug-specific: {drug_spec}")
CLASSIFY_PY

# ============================================================================
# Step 9: MultiQC report — CONVERGENCE #4 (QC + analysis)
# ============================================================================
echo "[Step 9] Running MultiQC (convergence #4)..."
if [ ! -f "${OUT}/multiqc/multiqc_report.html" ]; then
    multiqc \
        "${OUT}/fastqc_raw" "${OUT}/fastqc_trimmed" "${OUT}/trimmed" \
        -o "${OUT}/multiqc" \
        --force \
        2>&1 | tail -3
fi

# ============================================================================
# Step 10: Generate final report.csv
# ============================================================================
echo "[Step 10] Generating final report..."
python3 << 'REPORT_PY'
import csv
import os

OUT = os.environ.get("OUT", "outputs")
RESULTS = os.environ.get("RESULTS", "results")

report = []

# --- Count QC metrics ---
with open(f"{OUT}/comparison/count_qc.tsv") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            report.append((parts[0], parts[1]))

# --- Read mageck count log for mapping stats ---
count_log = f"{OUT}/count/screen.count_normalized.txt"
if os.path.exists(f"{OUT}/count/screen.countsummary.txt"):
    with open(f"{OUT}/count/screen.countsummary.txt") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            label = row.get('Label', '')
            reads = row.get('Reads', '0')
            mapped = row.get('Mapped', '0')
            pct = row.get('Percentage', '0')
            report.append((f"mapped_reads_{label.lower()}", mapped))
            report.append((f"mapping_pct_{label.lower()}", pct))

# --- RRA Drug results ---
with open(f"{OUT}/rra_drug/drug_vs_t0.gene_summary.txt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    genes = list(reader)

# Top depleted gene (negative selection)
genes_neg = sorted(genes, key=lambda x: float(x['neg|rank']))
if genes_neg:
    report.append(("top_depleted_gene_rra", genes_neg[0]['id']))
    report.append(("top_depleted_fdr_rra", genes_neg[0]['neg|fdr']))
    report.append(("top_depleted_lfc_rra", genes_neg[0]['neg|lfc']))

# Top enriched gene (positive selection)
genes_pos = sorted(genes, key=lambda x: float(x['pos|rank']))
if genes_pos:
    report.append(("top_enriched_gene_rra", genes_pos[0]['id']))
    report.append(("top_enriched_fdr_rra", genes_pos[0]['pos|fdr']))

# Count significant genes
num_dep = sum(1 for g in genes if float(g['neg|fdr']) < 0.25)
num_enr = sum(1 for g in genes if float(g['pos|fdr']) < 0.25)
report.append(("num_depleted_genes_fdr25_rra", str(num_dep)))
report.append(("num_enriched_genes_fdr25_rra", str(num_enr)))

# --- MLE results ---
mle_file = f"{OUT}/mle/screen_mle.gene_summary.txt"
if os.path.exists(mle_file):
    with open(mle_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        mle_genes = list(reader)
    if mle_genes:
        # Sort by drug|fdr
        try:
            mle_sorted = sorted(mle_genes, key=lambda x: float(x.get('drug|fdr', 1)))
            report.append(("top_gene_mle", mle_sorted[0].get('Gene', 'NA')))
            report.append(("top_gene_mle_fdr", mle_sorted[0].get('drug|fdr', 'NA')))
            report.append(("top_gene_mle_beta", mle_sorted[0].get('drug|beta', 'NA')))
            num_mle_sig = sum(1 for g in mle_genes if float(g.get('drug|fdr', 1)) < 0.25)
            report.append(("num_significant_mle_fdr25", str(num_mle_sig)))
        except (ValueError, KeyError):
            pass

# --- Vehicle results ---
with open(f"{OUT}/rra_vehicle/vehicle_vs_t0.gene_summary.txt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    veh_genes = list(reader)
num_veh_dep = sum(1 for g in veh_genes if float(g['neg|fdr']) < 0.25)
report.append(("num_depleted_genes_vehicle_fdr25", str(num_veh_dep)))

# --- Concordance ---
with open(f"{OUT}/comparison/classified_hits.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    classified = list(reader)
high_conf = sum(1 for h in classified if h['tier'] == 'high_confidence')
moderate_conf = sum(1 for h in classified if h['tier'] == 'moderate_confidence')
drug_specific = sum(1 for h in classified if h['drug_specific'] == 'True')
report.append(("high_confidence_hits", str(high_conf)))
report.append(("moderate_confidence_hits", str(moderate_conf)))
report.append(("drug_specific_depleted_genes", str(drug_specific)))

# --- Drug-specific top gene ---
with open(f"{OUT}/comparison/drug_specific_hits.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    drug_hits = list(reader)
if drug_hits:
    report.append(("top_drug_specific_gene", drug_hits[0]['gene']))
    report.append(("top_drug_specific_fdr", drug_hits[0]['rra_neg_fdr']))

# Write final report
with open(f"{RESULTS}/report.csv", 'w') as f:
    f.write("metric,value\n")
    for m, v in report:
        f.write(f"{m},{v}\n")

print(f"Report written with {len(report)} metrics")
for m, v in report:
    print(f"  {m}: {v}")
REPORT_PY

echo ""
echo "========================================="
echo "  CRISPR Screen Analysis Complete!"
echo "========================================="
echo ""
cat "${RESULTS}/report.csv"
