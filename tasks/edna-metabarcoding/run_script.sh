#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# eDNA Metabarcoding Pipeline: Aquatic Biodiversity Assessment
# ============================================================
# DAG structure (depth=10, convergence=4):
#
#  sample_R1.fq.gz + sample_R2.fq.gz (x6 samples)
#         │
#  [cutadapt primer removal] ─── per sample         Level 1
#         │
#  ┌──────┼──────────┐
#  │      │          │
# [vsearch  [fastqc  [seqkit                        Level 2
#  merge     QC]      stats]
#  pairs]
#  │      │          │
#  └──────┼──────────┘
#         │
#  [CONVERGENCE 1: QC + merged reads]               Level 3
#         │
#  [vsearch quality filter]                          Level 4
#         │
#  [pool samples + vsearch dereplicate]              Level 5
#         │
#  ┌──────┼──────────┐
#  │      │          │
# [vsearch [swarm   [vsearch                        Level 6
#  cluster  cluster] denoise
#  97%]              UNOISE3]
#  │      │          │
#  └──────┼──────────┘
#         │
#  [CONVERGENCE 2: select consensus + chimera removal] Level 7
#         │
#  ┌──────┼──────────┐
#  │                 │
# [BLAST           [vsearch                          Level 8
#  taxonomy]        usearch_global
#                   taxonomy]
#  │                 │
#  └────────┬────────┘
#           │
#  [CONVERGENCE 3: LCA consensus taxonomy]           Level 9
#           │
#  ┌────────┼──────────┐
#  │        │          │
# [species [R/vegan   [detection                     Level 9
#  list]    diversity]  probability]
#  │        │          │
#  └────────┼──────────┘
#           │
#  [CONVERGENCE 4: final report with QC]             Level 10
#
# Longest path: primer_removal -> merge -> QC_convergence ->
#   quality_filter -> dereplicate -> cluster -> chimera_removal ->
#   BLAST -> LCA -> diversity -> report = depth 10
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{trimmed,merged,filtered,qc,derep,clusters,chimera,taxonomy,community}
mkdir -p "${RESULTS}"

# Sample list
SAMPLES=(DRR205394 DRR205395 DRR205396 DRR205397 DRR205398 DRR205399)

# MiFish-U primer sequences (Miya et al. 2015)
FWD_PRIMER="GTCGGTAAAACTCGTGCCAGC"
REV_PRIMER="CATAGTGGGGTATCTAATCCCAGTTTG"
# Reverse complement of reverse primer
REV_PRIMER_RC=$(echo "$REV_PRIMER" | tr ACGTacgt TGCAtgca | rev)

# ============================================================
# Level 1: Primer removal with cutadapt (per sample)
# ============================================================
echo "=== Level 1: Primer removal ==="
for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/trimmed/${S}_R1.fastq.gz" ]; then
    cutadapt \
      -g "${FWD_PRIMER}" \
      -G "${REV_PRIMER}" \
      --discard-untrimmed \
      --minimum-length 50 \
      -j ${THREADS} \
      -o "${OUT}/trimmed/${S}_R1.fastq.gz" \
      -p "${OUT}/trimmed/${S}_R2.fastq.gz" \
      "${DATA}/${S}_R1.fastq.gz" \
      "${DATA}/${S}_R2.fastq.gz" \
      > "${OUT}/trimmed/${S}_cutadapt.log" 2>&1
    echo "  ${S}: trimmed"
  fi
done

# ============================================================
# Level 2: QC + merge + stats (parallel branches)
# ============================================================
echo "=== Level 2: QC + merge + stats ==="

# Branch 2a: FastQC on trimmed reads
if [ ! -f "${OUT}/qc/fastqc_done" ]; then
  for S in "${SAMPLES[@]}"; do
    fastqc -t ${THREADS} -o "${OUT}/qc/" \
      "${OUT}/trimmed/${S}_R1.fastq.gz" \
      "${OUT}/trimmed/${S}_R2.fastq.gz" \
      > /dev/null 2>&1
  done
  touch "${OUT}/qc/fastqc_done"
  echo "  FastQC done"
fi

# Branch 2b: seqkit stats on trimmed reads
if [ ! -f "${OUT}/qc/seqkit_stats.tsv" ]; then
  seqkit stats -T -j ${THREADS} "${OUT}/trimmed/"*.fastq.gz > "${OUT}/qc/seqkit_stats.tsv" 2>/dev/null
  echo "  seqkit stats done"
fi

# Branch 2c: vsearch merge pairs (per sample)
for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/merged/${S}.fastq" ]; then
    vsearch --fastq_mergepairs "${OUT}/trimmed/${S}_R1.fastq.gz" \
      --reverse "${OUT}/trimmed/${S}_R2.fastq.gz" \
      --fastqout "${OUT}/merged/${S}.fastq" \
      --fastq_maxdiffs 10 \
      --fastq_minovlen 50 \
      --threads ${THREADS} \
      --label_suffix ";sample=${S}" \
      > "${OUT}/merged/${S}_merge.log" 2>&1
    echo "  ${S}: merged"
  fi
done

# ============================================================
# Level 3: CONVERGENCE 1 — QC + merged reads available
# ============================================================
echo "=== Level 3: Convergence 1 (QC + merged) ==="
# MultiQC aggregation
if [ ! -f "${OUT}/qc/multiqc_report.html" ]; then
  multiqc "${OUT}/qc/" "${OUT}/trimmed/" -o "${OUT}/qc/" --force > /dev/null 2>&1 || true
  echo "  MultiQC done"
fi

# ============================================================
# Level 4: Quality filter (per sample)
# ============================================================
echo "=== Level 4: Quality filtering ==="
for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/filtered/${S}.fasta" ]; then
    vsearch --fastq_filter "${OUT}/merged/${S}.fastq" \
      --fastq_maxee 1.0 \
      --fastq_minlen 100 \
      --fastq_maxlen 300 \
      --fastaout "${OUT}/filtered/${S}.fasta" \
      --relabel "${S}." \
      > "${OUT}/filtered/${S}_filter.log" 2>&1
    echo "  ${S}: filtered"
  fi
done

# ============================================================
# Level 5: Pool samples + dereplicate
# ============================================================
echo "=== Level 5: Pool + dereplicate ==="
if [ ! -f "${OUT}/derep/all_derep.fasta" ]; then
  # Pool all filtered sequences
  cat "${OUT}/filtered/"*.fasta > "${OUT}/derep/all_pooled.fasta"

  # Dereplicate
  vsearch --derep_fulllength "${OUT}/derep/all_pooled.fasta" \
    --output "${OUT}/derep/all_derep.fasta" \
    --sizein --sizeout \
    --minuniquesize 2 \
    --uc "${OUT}/derep/all_derep.uc" \
    > "${OUT}/derep/derep.log" 2>&1
  echo "  Dereplication done"
fi

UNIQUE_COUNT=$(grep -c "^>" "${OUT}/derep/all_derep.fasta" || true)
echo "  Unique sequences: ${UNIQUE_COUNT}"

# ============================================================
# Level 6: Three parallel clustering methods
# ============================================================
echo "=== Level 6: Clustering (3 methods) ==="

# Method 6a: vsearch OTU clustering at 97%
if [ ! -f "${OUT}/clusters/otu97_centroids.fasta" ]; then
  vsearch --cluster_size "${OUT}/derep/all_derep.fasta" \
    --id 0.97 \
    --centroids "${OUT}/clusters/otu97_centroids.fasta" \
    --uc "${OUT}/clusters/otu97.uc" \
    --sizein --sizeout \
    --threads ${THREADS} \
    > "${OUT}/clusters/otu97.log" 2>&1
  echo "  OTU 97% clustering done"
fi

# Method 6b: SWARM clustering
if [ ! -f "${OUT}/clusters/swarm_centroids.fasta" ]; then
  # swarm needs dereplicated sequences sorted by abundance
  vsearch --sortbysize "${OUT}/derep/all_derep.fasta" \
    --output "${OUT}/clusters/sorted_for_swarm.fasta" \
    --sizein --sizeout 2>/dev/null

  swarm -d 1 -z \
    -w "${OUT}/clusters/swarm_centroids.fasta" \
    -o "${OUT}/clusters/swarm_otus.txt" \
    -s "${OUT}/clusters/swarm_stats.txt" \
    -t ${THREADS} \
    "${OUT}/clusters/sorted_for_swarm.fasta" \
    > "${OUT}/clusters/swarm.log" 2>&1
  echo "  SWARM clustering done"
fi

# Method 6c: vsearch UNOISE3 denoising (ASVs)
if [ ! -f "${OUT}/clusters/unoise3_asvs.fasta" ]; then
  vsearch --cluster_unoise "${OUT}/derep/all_derep.fasta" \
    --centroids "${OUT}/clusters/unoise3_asvs.fasta" \
    --sizein --sizeout \
    --minsize 2 \
    > "${OUT}/clusters/unoise3.log" 2>&1
  echo "  UNOISE3 denoising done"
fi

OTU97_COUNT=$(grep -c "^>" "${OUT}/clusters/otu97_centroids.fasta" || true)
SWARM_COUNT=$(grep -c "^>" "${OUT}/clusters/swarm_centroids.fasta" || true)
UNOISE3_COUNT=$(grep -c "^>" "${OUT}/clusters/unoise3_asvs.fasta" || true)
echo "  OTU97: ${OTU97_COUNT}, SWARM: ${SWARM_COUNT}, UNOISE3: ${UNOISE3_COUNT}"

# ============================================================
# Level 7: CONVERGENCE 2 — Select consensus + chimera removal
# ============================================================
echo "=== Level 7: Convergence 2 (consensus + chimera removal) ==="

# Use UNOISE3 ASVs as primary (most conservative denoising method)
# Then apply chimera removal
if [ ! -f "${OUT}/chimera/clean_asvs.fasta" ]; then
  vsearch --uchime_denovo "${OUT}/clusters/unoise3_asvs.fasta" \
    --nonchimeras "${OUT}/chimera/clean_asvs.fasta" \
    --chimeras "${OUT}/chimera/chimeras.fasta" \
    --sizein --sizeout \
    > "${OUT}/chimera/chimera.log" 2>&1
  echo "  Chimera removal done"
fi

CLEAN_COUNT=$(grep -c "^>" "${OUT}/chimera/clean_asvs.fasta" || true)
CHIMERA_COUNT=$(grep -c "^>" "${OUT}/chimera/chimeras.fasta" 2>/dev/null || true)
CHIMERA_COUNT=${CHIMERA_COUNT:-0}
echo "  Clean ASVs: ${CLEAN_COUNT}, Chimeras removed: ${CHIMERA_COUNT}"

# ============================================================
# Level 7.5: Build BLAST database from reference
# ============================================================
echo "=== Building BLAST database ==="
if [ ! -f "${REF}/mitofish_12S.ndb" ]; then
  makeblastdb -in "${REF}/mitofish_12S.fasta" \
    -dbtype nucl \
    -out "${REF}/mitofish_12S" \
    -parse_seqids \
    > /dev/null 2>&1
  echo "  BLAST DB built"
fi

# ============================================================
# Level 8: Taxonomy assignment (two parallel methods)
# ============================================================
echo "=== Level 8: Taxonomy assignment ==="

# Method 8a: BLAST against MitoFish 12S
if [ ! -f "${OUT}/taxonomy/blast_hits.tsv" ]; then
  # Strip size annotations from headers for BLAST
  sed 's/;size=[0-9]*//' "${OUT}/chimera/clean_asvs.fasta" > "${OUT}/taxonomy/query.fasta"

  blastn -query "${OUT}/taxonomy/query.fasta" \
    -db "${REF}/mitofish_12S" \
    -out "${OUT}/taxonomy/blast_hits.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 1e-10 \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    > /dev/null 2>&1
  echo "  BLAST done"
fi

# Method 8b: vsearch usearch_global against MitoFish 12S
if [ ! -f "${OUT}/taxonomy/vsearch_hits.tsv" ]; then
  vsearch --usearch_global "${OUT}/taxonomy/query.fasta" \
    --db "${REF}/mitofish_12S.fasta" \
    --id 0.80 \
    --maxaccepts 10 \
    --blast6out "${OUT}/taxonomy/vsearch_hits.tsv" \
    --threads ${THREADS} \
    > "${OUT}/taxonomy/vsearch.log" 2>&1
  echo "  vsearch global search done"
fi

# ============================================================
# Level 9: CONVERGENCE 3 — LCA consensus taxonomy
# ============================================================
echo "=== Level 9: Convergence 3 (LCA taxonomy) ==="

if [ ! -f "${OUT}/taxonomy/lca_taxonomy.tsv" ]; then
  python3 << 'PYEOF'
import csv
import sys
from collections import defaultdict

# Load MitoFish taxonomy lookup
tax_lookup = {}
with open("reference/mitofish_12S_taxonomy.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        acc = row['Accession']
        tax_lookup[acc] = {
            'superkingdom': row.get('Superkingdom', ''),
            'phylum': row.get('Phylum', ''),
            'class': row.get('Class', ''),
            'order': row.get('Order', ''),
            'family': row.get('Family', ''),
            'genus': row.get('Genus', ''),
            'species': row.get('Species', '')
        }

# Read BLAST hits (top hits per query)
blast_tax = defaultdict(list)
with open("outputs/taxonomy/blast_hits.tsv") as f:
    for line in f:
        parts = line.strip().split('\t')
        qid, sid, pident = parts[0], parts[1], float(parts[2])
        if pident >= 97.0 and sid in tax_lookup:
            blast_tax[qid].append(tax_lookup[sid])

# Read vsearch hits
vsearch_tax = defaultdict(list)
with open("outputs/taxonomy/vsearch_hits.tsv") as f:
    for line in f:
        parts = line.strip().split('\t')
        qid, sid, pident = parts[0], parts[1], float(parts[2])
        if pident >= 97.0 and sid in tax_lookup:
            vsearch_tax[qid].append(tax_lookup[sid])

# LCA function
RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
def lca(tax_list):
    if not tax_list:
        return {r: '' for r in RANKS}
    result = {}
    for rank in RANKS:
        values = set(t[rank] for t in tax_list if t[rank])
        if len(values) == 1:
            result[rank] = values.pop()
        else:
            result[rank] = ''
            break  # stop at first disagreement
    for rank in RANKS:
        if rank not in result:
            result[rank] = ''
    return result

# Merge BLAST + vsearch via LCA
all_queries = set(list(blast_tax.keys()) + list(vsearch_tax.keys()))
with open("outputs/taxonomy/lca_taxonomy.tsv", 'w') as f:
    f.write("asv_id\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for qid in sorted(all_queries):
        combined = blast_tax.get(qid, []) + vsearch_tax.get(qid, [])
        tax = lca(combined)
        f.write(f"{qid}\t{tax['superkingdom']}\t{tax['phylum']}\t{tax['class']}\t{tax['order']}\t{tax['family']}\t{tax['genus']}\t{tax['species']}\n")

print(f"LCA taxonomy assigned to {len(all_queries)} ASVs")
PYEOF
fi

# ============================================================
# Level 9 continued: Three parallel community analyses
# ============================================================
echo "=== Level 9: Community analyses ==="

# Build OTU table (ASV x sample) by mapping reads back
if [ ! -f "${OUT}/community/otu_table.tsv" ]; then
  # Map all filtered reads back to clean ASVs
  vsearch --usearch_global "${OUT}/derep/all_pooled.fasta" \
    --db "${OUT}/chimera/clean_asvs.fasta" \
    --id 0.97 \
    --otutabout "${OUT}/community/otu_table.tsv" \
    --threads ${THREADS} \
    > "${OUT}/community/map.log" 2>&1
  echo "  OTU table built"
fi

# Branch 9a: Species list
# Branch 9b: Alpha + beta diversity (R/vegan)
# Branch 9c: Detection probability
if [ ! -f "${OUT}/community/diversity_results.tsv" ]; then
  python3 << 'PYEOF'
import csv
import math
from collections import defaultdict

# Load OTU table
otu_table = {}  # {asv_id: {sample: count}}
samples = []
with open("outputs/community/otu_table.tsv") as f:
    header = f.readline().strip().split('\t')
    samples = header[1:]  # first col is OTU ID
    for line in f:
        parts = line.strip().split('\t')
        asv_id = parts[0]
        counts = [int(x) for x in parts[1:]]
        otu_table[asv_id] = dict(zip(samples, counts))

# Load taxonomy
taxonomy = {}
with open("outputs/taxonomy/lca_taxonomy.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        taxonomy[row['asv_id']] = row

# === Branch 9a: Species list ===
species_set = set()
genus_set = set()
family_set = set()
order_set = set()
for asv_id, tax in taxonomy.items():
    if tax['species']:
        species_set.add(tax['species'])
    if tax['genus']:
        genus_set.add(tax['genus'])
    if tax['family']:
        family_set.add(tax['family'])
    if tax['order']:
        order_set.add(tax['order'])

with open("outputs/community/species_list.tsv", 'w') as f:
    f.write("species\tgenus\tfamily\torder\n")
    for sp in sorted(species_set):
        # find matching taxonomy
        for asv_id, tax in taxonomy.items():
            if tax['species'] == sp:
                f.write(f"{sp}\t{tax['genus']}\t{tax['family']}\t{tax['order']}\n")
                break

# === Branch 9b: Alpha diversity (Shannon index per sample) ===
shannon_per_sample = {}
richness_per_sample = {}
for s in samples:
    counts = [otu_table[asv][s] for asv in otu_table if otu_table[asv].get(s, 0) > 0]
    total = sum(counts)
    if total == 0:
        shannon_per_sample[s] = 0.0
        richness_per_sample[s] = 0
        continue
    richness_per_sample[s] = len(counts)
    shannon = 0.0
    for c in counts:
        p = c / total
        if p > 0:
            shannon -= p * math.log(p)
    shannon_per_sample[s] = round(shannon, 4)

with open("outputs/community/diversity_results.tsv", 'w') as f:
    f.write("sample\tshannon_diversity\tspecies_richness\n")
    for s in samples:
        f.write(f"{s}\t{shannon_per_sample[s]}\t{richness_per_sample[s]}\n")

# === Branch 9c: Detection probability ===
# For each species, proportion of samples where detected
species_detection = defaultdict(int)
species_total_reads = defaultdict(int)
for asv_id in otu_table:
    sp = taxonomy.get(asv_id, {}).get('species', '')
    if not sp:
        continue
    for s in samples:
        if otu_table[asv_id].get(s, 0) > 0:
            species_detection[sp] += 1
            species_total_reads[sp] += otu_table[asv_id][s]

with open("outputs/community/detection_probability.tsv", 'w') as f:
    f.write("species\tsamples_detected\tdetection_rate\ttotal_reads\n")
    for sp in sorted(species_detection.keys()):
        det_rate = round(species_detection[sp] / len(samples), 4)
        f.write(f"{sp}\t{species_detection[sp]}\t{det_rate}\t{species_total_reads[sp]}\n")

print(f"Species: {len(species_set)}, Genera: {len(genus_set)}, Families: {len(family_set)}, Orders: {len(order_set)}")
print(f"Shannon range: {min(shannon_per_sample.values()):.4f} - {max(shannon_per_sample.values()):.4f}")
PYEOF
  echo "  Python community analysis done"
fi

# R/vegan beta diversity
if [ ! -f "${OUT}/community/beta_diversity.tsv" ]; then
  cat > "${OUT}/community/run_vegan.R" << 'REOF'
library(vegan)
otu <- read.delim("outputs/community/otu_table.tsv", row.names=1, check.names=FALSE)
otu_t <- t(otu)
bc <- as.matrix(vegdist(otu_t, method="bray"))
write.table(bc, "outputs/community/beta_diversity.tsv", sep="\t", quote=FALSE)
cat("Beta diversity (Bray-Curtis) range:", range(bc[lower.tri(bc)]), "\n")
cat("Mean Bray-Curtis:", mean(bc[lower.tri(bc)]), "\n")
REOF
  Rscript "${OUT}/community/run_vegan.R"
  echo "  R/vegan beta diversity done"
fi

# ============================================================
# Level 10: CONVERGENCE 4 — Final report
# ============================================================
echo "=== Level 10: Final report ==="

python3 << 'PYEOF'
import csv
import math
import os
from collections import defaultdict

# Gather all metrics
metrics = {}

# --- Raw read counts ---
total_raw = 0
for s in ["DRR205394","DRR205395","DRR205396","DRR205397","DRR205398","DRR205399"]:
    for r in ["R1", "R2"]:
        fpath = f"outputs/qc/seqkit_stats.tsv"
        break
    break

# Count from seqkit stats
with open("outputs/qc/seqkit_stats.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        total_raw += int(row['num_seqs'].replace(',', ''))
metrics['total_raw_reads'] = total_raw

# --- Merged read counts ---
total_merged = 0
for s in ["DRR205394","DRR205395","DRR205396","DRR205397","DRR205398","DRR205399"]:
    logf = f"outputs/merged/{s}_merge.log"
    if os.path.exists(logf):
        with open(logf) as f:
            for line in f:
                if "Merged" in line and "pairs" in line:
                    # Parse vsearch merge log
                    parts = line.strip().split()
                    for i, p in enumerate(parts):
                        if p.isdigit() or p.replace(',','').isdigit():
                            total_merged += int(p.replace(',',''))
                            break
                    break

# Fallback: count merged reads directly
if total_merged == 0:
    for s in ["DRR205394","DRR205395","DRR205396","DRR205397","DRR205398","DRR205399"]:
        mf = f"outputs/merged/{s}.fastq"
        if os.path.exists(mf):
            count = sum(1 for line in open(mf)) // 4
            total_merged += count
metrics['total_merged_reads'] = total_merged

# Merge rate
if total_raw > 0:
    # total_raw is R1+R2, so pairs = total_raw / 2
    metrics['merge_rate'] = round(total_merged / (total_raw / 2) * 100, 2)
else:
    metrics['merge_rate'] = 0.0

# --- Unique sequences ---
derep_fasta = "outputs/derep/all_derep.fasta"
unique_count = sum(1 for line in open(derep_fasta) if line.startswith('>'))
metrics['unique_sequences'] = unique_count

# --- Clustering results ---
for method, fname in [("clusters_otu97", "outputs/clusters/otu97_centroids.fasta"),
                       ("clusters_swarm", "outputs/clusters/swarm_centroids.fasta"),
                       ("clusters_denoised", "outputs/clusters/unoise3_asvs.fasta")]:
    count = sum(1 for line in open(fname) if line.startswith('>'))
    metrics[method] = count

# --- Chimera removal ---
clean_fasta = "outputs/chimera/clean_asvs.fasta"
chimera_fasta = "outputs/chimera/chimeras.fasta"
metrics['clean_sequence_count'] = sum(1 for line in open(clean_fasta) if line.startswith('>'))
chimera_count = 0
if os.path.exists(chimera_fasta):
    chimera_count = sum(1 for line in open(chimera_fasta) if line.startswith('>'))
metrics['chimera_count'] = chimera_count

# --- Taxonomy assignment ---
assigned = 0
total_asvs = 0
with open("outputs/taxonomy/lca_taxonomy.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        total_asvs += 1
        if row['species'] or row['genus'] or row['family']:
            assigned += 1
metrics['assigned_sequences'] = assigned
metrics['unassigned_sequences'] = total_asvs - assigned

# --- Species/genus/family/order counts ---
species_set = set()
genus_set = set()
family_set = set()
order_set = set()
with open("outputs/taxonomy/lca_taxonomy.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['species']: species_set.add(row['species'])
        if row['genus']: genus_set.add(row['genus'])
        if row['family']: family_set.add(row['family'])
        if row['order']: order_set.add(row['order'])
metrics['species_count'] = len(species_set)
metrics['genus_count'] = len(genus_set)
metrics['family_count'] = len(family_set)
metrics['order_count'] = len(order_set)

# --- Diversity ---
shannon_vals = []
richness_vals = []
with open("outputs/community/diversity_results.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        shannon_vals.append(float(row['shannon_diversity']))
        richness_vals.append(int(row['species_richness']))
metrics['mean_shannon_diversity'] = round(sum(shannon_vals)/len(shannon_vals), 4)
metrics['min_species_richness'] = min(richness_vals)
metrics['max_species_richness'] = max(richness_vals)

# --- Beta diversity ---
bc_vals = []
with open("outputs/community/beta_diversity.tsv") as f:
    header = f.readline().strip().split('\t')
    rows = []
    for line in f:
        parts = line.strip().split('\t')
        rows.append([float(x) for x in parts[1:]])
    for i in range(len(rows)):
        for j in range(i+1, len(rows)):
            bc_vals.append(rows[i][j])
if bc_vals:
    metrics['mean_beta_diversity'] = round(sum(bc_vals)/len(bc_vals), 4)
    metrics['min_beta_diversity'] = round(min(bc_vals), 4)
    metrics['max_beta_diversity'] = round(max(bc_vals), 4)

# --- Detection rate ---
det_rates = []
with open("outputs/community/detection_probability.tsv") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        det_rates.append(float(row['detection_rate']))
if det_rates:
    metrics['mean_detection_rate'] = round(sum(det_rates)/len(det_rates), 4)

# === Write report ===
with open("results/report.csv", 'w') as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("=== Report generated ===")
for k, v in metrics.items():
    print(f"  {k} = {v}")
PYEOF

echo "=== Pipeline complete ==="
