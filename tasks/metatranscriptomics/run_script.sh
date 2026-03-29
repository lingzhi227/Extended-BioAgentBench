#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Metatranscriptomics Pipeline — DAG (depth=10, convergence=4)
# ============================================================
#
#  R1.fastq.gz     R2.fastq.gz
#      │                │
#  [fastp QC] ────── [fastp QC]                     Level 1
#      │                │
#      └───────┬────────┘
#              │
#  [bbduk rRNA removal] ◄── rRNA kmer databases     Level 2
#              │
#      ┌───────┴──────────────┐
#      │                      │
#  (non-rRNA reads)       (rRNA reads)
#      │                      │
#      │              [python rRNA                   Level 3
#      │               community profile]
#      │                      │
#  ┌───┴────────────┐         │
#  │                │         │
# [MEGAHIT      [samtools     │                      Level 4
#  assembly]     flagstat     │
#  │             (read stats)]│
#  │                │         │
# [prodigal         │         │                      Level 5
#  gene calling]    │         │
#  │       │        │         │
# [diamond [bowtie2 │         │                      Level 6
#  blastx]  map to  │         │
#           contigs]│         │
#  │       │        │         │
#  │   [samtools    │         │
#  │    sort+index] │         │
#  │       │        │         │
#  │   [featureCounts]        │                      Level 7
#  │       │        │         │
#  └───────┴────────┘         │
#              │              │
#      [CONVERGENCE 1] ◄─────┘                      Level 7
#      (assembly+func+counts+rRNA taxonomy)
#              │
#      ┌───────┼───────────┐
#      │       │           │
#  [python  [python     [python                      Level 8
#   func     expression   taxonomy
#   summary] analysis]    merge]
#      │       │           │
#      └───────┼───────────┘
#              │
#      [CONVERGENCE 2]                               Level 8
#      (func groups + expression + taxonomy)
#              │
#      ┌───────┼───────────┐
#      │       │           │
#  [python  [python     [python                      Level 9
#   diversity active      functional
#   metrics]  taxa]       enrichment]
#      │       │           │
#      └───────┼───────────┘
#              │
#      [CONVERGENCE 3]                               Level 9
#      (diversity + active taxa + enrichment)
#              │
#      [CONVERGENCE 4] ◄── QC stats + rRNA%          Level 10
#      [python report]
#
# Convergence points:
#   C1: assembly + functional annotation + gene counts + rRNA taxonomy
#   C2: functional summary + expression + taxonomy
#   C3: diversity + active taxa + enrichment
#   C4: final report combining all + QC metrics
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORK=$(pwd)
DATA="${WORK}/data"
REF="${WORK}/reference"
OUT="${WORK}/outputs"
RESULTS="${WORK}/results"

mkdir -p "${OUT}"/{qc,rrna_filter,assembly,mapping,genes,annotation,counts,analysis} "${RESULTS}"

# ─── Level 1: Read QC with fastp ───
if [ ! -f "${OUT}/qc/reads_R1.trimmed.fastq.gz" ]; then
  echo "[Level 1] Running fastp QC..."
  fastp \
    -i "${DATA}/reads_R1.fastq.gz" \
    -I "${DATA}/reads_R2.fastq.gz" \
    -o "${OUT}/qc/reads_R1.trimmed.fastq.gz" \
    -O "${OUT}/qc/reads_R2.trimmed.fastq.gz" \
    --json "${OUT}/qc/fastp.json" \
    --html "${OUT}/qc/fastp.html" \
    --thread ${THREADS} \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --detect_adapter_for_pe
fi

# Extract QC stats
READS_BEFORE=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
READS_AFTER=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
Q30_BEFORE=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(round(d['summary']['before_filtering']['q30_rate']*100,2))")
Q30_AFTER=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

echo "  Reads before: ${READS_BEFORE}, after: ${READS_AFTER}"

# ─── Level 2: rRNA removal via alignment to rRNA databases ───
if [ ! -f "${OUT}/rrna_filter/non_rrna_R1.fastq.gz" ]; then
  echo "[Level 2] Running rRNA removal..."

  # Concatenate rRNA reference databases
  cat "${REF}"/rrna_db/*.fasta > "${OUT}/rrna_filter/rrna_combined.fasta"

  # Build bowtie2 index from rRNA references
  bowtie2-build \
    "${OUT}/rrna_filter/rrna_combined.fasta" \
    "${OUT}/rrna_filter/rrna_index" \
    --threads ${THREADS} --quiet

  # Align reads to rRNA — extract mapped (rRNA) reads
  bowtie2 \
    -x "${OUT}/rrna_filter/rrna_index" \
    -1 "${OUT}/qc/reads_R1.trimmed.fastq.gz" \
    -2 "${OUT}/qc/reads_R2.trimmed.fastq.gz" \
    --threads ${THREADS} \
    --very-fast \
    --no-unal \
    2> "${OUT}/rrna_filter/bowtie2_rrna.log" \
  | samtools sort -n -@ 4 -o "${OUT}/rrna_filter/rrna_aligned.bam"

  # Extract rRNA reads as FASTQ
  samtools fastq \
    -1 "${OUT}/rrna_filter/rrna_R1.fastq.gz" \
    -2 "${OUT}/rrna_filter/rrna_R2.fastq.gz" \
    -0 /dev/null -s /dev/null \
    -F 4 \
    "${OUT}/rrna_filter/rrna_aligned.bam"

  # Extract non-rRNA (unmapped) reads
  bowtie2 \
    -x "${OUT}/rrna_filter/rrna_index" \
    -1 "${OUT}/qc/reads_R1.trimmed.fastq.gz" \
    -2 "${OUT}/qc/reads_R2.trimmed.fastq.gz" \
    --threads ${THREADS} \
    --very-fast \
    --un-conc-gz "${OUT}/rrna_filter/non_rrna_%.fastq.gz" \
    2> /dev/null \
  | samtools view -c > /dev/null

  mv "${OUT}/rrna_filter/non_rrna_1.fastq.gz" "${OUT}/rrna_filter/non_rrna_R1.fastq.gz"
  mv "${OUT}/rrna_filter/non_rrna_2.fastq.gz" "${OUT}/rrna_filter/non_rrna_R2.fastq.gz"
fi

# Count rRNA vs non-rRNA reads
RRNA_READS_R1=$(zcat "${OUT}/rrna_filter/rrna_R1.fastq.gz" 2>/dev/null | awk 'NR%4==1' | wc -l || true)
RRNA_READS=$(( RRNA_READS_R1 * 2 ))
NONRRNA_READS_R1=$(zcat "${OUT}/rrna_filter/non_rrna_R1.fastq.gz" 2>/dev/null | awk 'NR%4==1' | wc -l || true)
NONRRNA_READS=$(( NONRRNA_READS_R1 * 2 ))
TOTAL_FILTERED=$((RRNA_READS + NONRRNA_READS))
RRNA_PCT=$(python3 -c "print(round(${RRNA_READS}/${TOTAL_FILTERED}*100,2) if ${TOTAL_FILTERED}>0 else 0)")
echo "  rRNA reads: ${RRNA_READS} (${RRNA_PCT}%), non-rRNA: ${NONRRNA_READS}"

# ─── Level 3: rRNA community profile (parallel with assembly path) ───
if [ ! -f "${OUT}/analysis/rrna_community.tsv" ]; then
  echo "[Level 3] Profiling rRNA community..."
  python3 << 'PYEOF'
import gzip, os

out = os.environ.get("OUT", "outputs")
os.makedirs(f"{out}/analysis", exist_ok=True)

# Analyze rRNA reads: GC content, length distribution
seq_count = 0
gc_sum = 0
len_sum = 0
len_dist = {}
with gzip.open(f"{out}/rrna_filter/rrna_R1.fastq.gz", "rt") as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            seq = line.strip()
            seq_count += 1
            gc_sum += seq.count("G") + seq.count("C")
            slen = len(seq)
            len_sum += slen
            bucket = (slen // 25) * 25
            len_dist[bucket] = len_dist.get(bucket, 0) + 1

avg_gc = round(gc_sum / len_sum * 100, 2) if len_sum > 0 else 0
avg_len = round(len_sum / seq_count, 1) if seq_count > 0 else 0

with open(f"{out}/analysis/rrna_community.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"rrna_read_count\t{seq_count}\n")
    f.write(f"rrna_avg_gc_pct\t{avg_gc}\n")
    f.write(f"rrna_avg_length\t{avg_len}\n")
    for bucket in sorted(len_dist):
        f.write(f"rrna_len_{bucket}_{bucket+24}\t{len_dist[bucket]}\n")

print(f"  rRNA community: {seq_count} reads, avg GC={avg_gc}%, avg len={avg_len}")
PYEOF
fi

# ─── Level 4: MEGAHIT de novo assembly of non-rRNA reads ───
if [ ! -f "${OUT}/assembly/final.contigs.fa" ]; then
  echo "[Level 4] Running MEGAHIT assembly..."
  rm -rf "${OUT}/assembly"  # MEGAHIT fails if dir exists (M6)
  megahit \
    -1 "${OUT}/rrna_filter/non_rrna_R1.fastq.gz" \
    -2 "${OUT}/rrna_filter/non_rrna_R2.fastq.gz" \
    -o "${OUT}/assembly" \
    --min-contig-len 200 \
    -t ${THREADS}
fi

# Assembly stats
TOTAL_CONTIGS=$(grep -c "^>" "${OUT}/assembly/final.contigs.fa" || true)
ASSEMBLY_LENGTH=$(awk '/^>/{next}{sum+=length($0)}END{print sum}' "${OUT}/assembly/final.contigs.fa")
LARGEST_CONTIG=$(awk '/^>/{if(len>max)max=len;len=0;next}{len+=length($0)}END{if(len>max)max=len;print max}' "${OUT}/assembly/final.contigs.fa")
ASSEMBLY_GC=$(awk '/^>/{next}{for(i=1;i<=length($0);i++){c=substr($0,i,1);if(c=="G"||c=="C"||c=="g"||c=="c")gc++;tot++}}END{printf "%.2f",gc/tot*100}' "${OUT}/assembly/final.contigs.fa")
N50=$(python3 -c "
lens=[]
with open('${OUT}/assembly/final.contigs.fa') as f:
    l=0
    for line in f:
        if line.startswith('>'):
            if l: lens.append(l)
            l=0
        else: l+=len(line.strip())
    if l: lens.append(l)
lens.sort(reverse=True)
total=sum(lens)
cum=0
for x in lens:
    cum+=x
    if cum>=total/2:
        print(x); break
")
echo "  Assembly: ${TOTAL_CONTIGS} contigs, ${ASSEMBLY_LENGTH} bp, N50=${N50}"

# ─── Level 5: Gene prediction with prodigal ───
if [ ! -f "${OUT}/genes/genes.gff" ]; then
  echo "[Level 5] Running prodigal gene prediction..."
  prodigal \
    -i "${OUT}/assembly/final.contigs.fa" \
    -o "${OUT}/genes/genes.gff" \
    -a "${OUT}/genes/proteins.faa" \
    -d "${OUT}/genes/genes.fna" \
    -p meta \
    -f gff
fi
PREDICTED_GENES=$(grep -c "^>" "${OUT}/genes/proteins.faa" || true)
echo "  Predicted genes: ${PREDICTED_GENES}"

# ─── Level 6a: DIAMOND blastx functional annotation ───
if [ ! -f "${OUT}/annotation/diamond_hits.tsv" ]; then
  echo "[Level 6a] Running DIAMOND blastx..."
  diamond blastx \
    --query "${OUT}/genes/genes.fna" \
    --db "${REF}/uniprot_sprot.dmnd" \
    --out "${OUT}/annotation/diamond_hits.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --threads ${THREADS} \
    --sensitive
fi
ANNOTATED_GENES=$(cut -f1 "${OUT}/annotation/diamond_hits.tsv" | sort -u | wc -l || true)
echo "  Annotated genes: ${ANNOTATED_GENES}"

# ─── Level 6b: Map reads back to assembled contigs with bowtie2 ───
if [ ! -f "${OUT}/mapping/mapped.sorted.bam" ]; then
  echo "[Level 6b] Building bowtie2 index and mapping..."
  bowtie2-build \
    "${OUT}/assembly/final.contigs.fa" \
    "${OUT}/mapping/contigs_index" \
    --threads ${THREADS} \
    --quiet

  bowtie2 \
    -x "${OUT}/mapping/contigs_index" \
    -1 "${OUT}/rrna_filter/non_rrna_R1.fastq.gz" \
    -2 "${OUT}/rrna_filter/non_rrna_R2.fastq.gz" \
    --no-unal \
    --threads ${THREADS} \
    --rg-id sample \
    --rg "SM:sample" \
    2> "${OUT}/mapping/bowtie2.log" \
  | samtools sort -@ ${THREADS} -o "${OUT}/mapping/mapped.sorted.bam"

  samtools index "${OUT}/mapping/mapped.sorted.bam"
fi

# Mapping stats via samtools flagstat
samtools flagstat "${OUT}/mapping/mapped.sorted.bam" > "${OUT}/mapping/flagstat.txt"
MAPPED_READS=$(grep "mapped (" "${OUT}/mapping/flagstat.txt" | head -1 | awk '{print $1}')
# Mapping rate relative to total non-rRNA reads (not just reads in BAM, since --no-unal)
MAPPING_PCT=$(python3 -c "print(round(${MAPPED_READS}/${NONRRNA_READS}*100 if ${NONRRNA_READS}>0 else 0, 2))")
echo "  Mapped reads: ${MAPPED_READS}/${NONRRNA_READS} (${MAPPING_PCT}%)"

# ─── Level 7: featureCounts gene quantification ── CONVERGENCE 1 ───
if [ ! -f "${OUT}/counts/gene_counts.tsv" ]; then
  echo "[Level 7 / CONVERGENCE 1] Running featureCounts..."
  # Convert prodigal GFF to SAF format for featureCounts
  python3 << 'PYEOF'
import os
out = os.environ.get("OUT", "outputs")
with open(f"{out}/genes/genes.gff") as fin, open(f"{out}/counts/genes.saf", "w") as fout:
    fout.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for line in fin:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        chrom = parts[0]
        start = parts[3]
        end = parts[4]
        strand = parts[6]
        attrs = parts[8]
        gene_id = None
        for attr in attrs.split(";"):
            if attr.startswith("ID="):
                gene_id = attr.split("=")[1]
                break
        if gene_id:
            fout.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
PYEOF

  featureCounts \
    -a "${OUT}/counts/genes.saf" \
    -F SAF \
    -o "${OUT}/counts/gene_counts.tsv" \
    -p --countReadPairs \
    -T ${THREADS} \
    "${OUT}/mapping/mapped.sorted.bam"
fi

EXPRESSED_GENES=$(awk 'NR>2 && $NF>0{c++}END{print c+0}' "${OUT}/counts/gene_counts.tsv")
echo "  Expressed genes (count>0): ${EXPRESSED_GENES}"

# ─── Level 8: CONVERGENCE 2 — Functional summary + Expression analysis + Taxonomy merge ───
echo "[Level 8 / CONVERGENCE 2] Running integrated analyses..."
python3 << 'PYEOF'
import os, collections

out = os.environ.get("OUT", "outputs")

# --- Functional summary from DIAMOND hits ---
func_counts = collections.Counter()
org_counts = collections.Counter()
with open(f"{out}/annotation/diamond_hits.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 13:
            title = parts[12]
            func = title.split(" OS=")[0] if " OS=" in title else title
            # Strip Swiss-Prot accession prefix (e.g., "sp|Q46508|HNDD_SOLFR ")
            if func.startswith("sp|") or func.startswith("tr|"):
                parts2 = func.split(" ", 1)
                func = parts2[1] if len(parts2) > 1 else func
            func_counts[func] += 1
            if " OS=" in title:
                org = title.split(" OS=")[1].split(" OX=")[0]
                org_counts[org] += 1

top_funcs = func_counts.most_common(20)
with open(f"{out}/analysis/functional_summary.tsv", "w") as f:
    f.write("function\tcount\n")
    for func, count in top_funcs:
        f.write(f"{func}\t{count}\n")

top_orgs = org_counts.most_common(20)
with open(f"{out}/analysis/organism_distribution.tsv", "w") as f:
    f.write("organism\thit_count\n")
    for org, count in top_orgs:
        f.write(f"{org}\t{count}\n")

# --- Expression analysis: CPM normalization ---
gene_counts = {}
total_count = 0
with open(f"{out}/counts/gene_counts.tsv") as f:
    for line in f:
        if line.startswith("#") or line.startswith("Geneid"):
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 7:
            gid = parts[0]
            count = int(parts[6])
            gene_counts[gid] = count
            total_count += count

with open(f"{out}/analysis/expression_cpm.tsv", "w") as f:
    f.write("gene_id\traw_count\tcpm\n")
    for gid, count in sorted(gene_counts.items(), key=lambda x: -x[1])[:100]:
        cpm = round(count / total_count * 1e6, 2) if total_count > 0 else 0
        f.write(f"{gid}\t{count}\t{cpm}\n")

print(f"  Functional: {len(func_counts)} unique functions, top={top_funcs[0][0][:50] if top_funcs else 'N/A'}")
print(f"  Organisms: {len(org_counts)} unique, top={top_orgs[0][0] if top_orgs else 'N/A'}")
print(f"  Expression: {sum(1 for v in gene_counts.values() if v > 0)} expressed, total counts={total_count}")
PYEOF

# ─── Level 9: CONVERGENCE 3 — Diversity + Active taxa + Functional enrichment ───
echo "[Level 9 / CONVERGENCE 3] Computing diversity and enrichment..."
python3 << 'PYEOF'
import os, math, collections

out = os.environ.get("OUT", "outputs")

# --- Shannon diversity from organism distribution ---
org_counts = {}
with open(f"{out}/analysis/organism_distribution.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            org_counts[parts[0]] = int(parts[1])

total = sum(org_counts.values())
shannon = 0
for count in org_counts.values():
    if count > 0 and total > 0:
        p = count / total
        shannon -= p * math.log(p)
shannon = round(shannon, 4)

richness = len(org_counts)

simpson = 0
for count in org_counts.values():
    if total > 1:
        simpson += (count * (count - 1)) / (total * (total - 1))
simpson = round(1 - simpson, 4)

singletons = sum(1 for c in org_counts.values() if c == 1)
doubletons = sum(1 for c in org_counts.values() if c == 2)
chao1 = richness + (singletons * (singletons - 1)) / (2 * (doubletons + 1))
chao1 = round(chao1, 1)

# Active taxa report
top_active = sorted(org_counts.items(), key=lambda x: -x[1])[:10]
with open(f"{out}/analysis/active_taxa.tsv", "w") as f:
    f.write("taxon\thit_count\trelative_abundance_pct\n")
    for org, count in top_active:
        pct = round(count / total * 100, 2) if total > 0 else 0
        f.write(f"{org}\t{count}\t{pct}\n")

# Functional enrichment (keyword frequency)
keyword_counts = collections.Counter()
with open(f"{out}/analysis/functional_summary.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            func = parts[0].lower()
            count = int(parts[1])
            for word in ["kinase", "transferase", "synthase", "reductase", "oxidase",
                        "dehydrogenase", "transporter", "permease", "lyase", "ligase",
                        "protease", "peptidase", "hydrolase", "isomerase", "helicase",
                        "polymerase", "ribosomal", "elongation", "translation",
                        "transcription", "membrane", "binding", "regulatory"]:
                if word in func:
                    keyword_counts[word] += count

with open(f"{out}/analysis/functional_enrichment.tsv", "w") as f:
    f.write("functional_category\tgene_count\n")
    for kw, count in keyword_counts.most_common(20):
        f.write(f"{kw}\t{count}\n")

# Save diversity
with open(f"{out}/analysis/diversity_metrics.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"shannon_diversity\t{shannon}\n")
    f.write(f"simpson_diversity\t{simpson}\n")
    f.write(f"observed_richness\t{richness}\n")
    f.write(f"chao1_estimate\t{chao1}\n")

print(f"  Shannon={shannon}, Simpson={simpson}, Richness={richness}, Chao1={chao1}")
print(f"  Top active: {top_active[0][0] if top_active else 'N/A'}")
print(f"  Functional categories: {len(keyword_counts)}")
PYEOF

# ─── Level 10: CONVERGENCE 4 — Final report ───
echo "[Level 10 / CONVERGENCE 4] Generating final report..."
python3 << PYEOF
import os

out = os.environ.get("OUT", "outputs")
results = os.environ.get("RESULTS", "results")
os.makedirs(results, exist_ok=True)

# Read diversity metrics
diversity = {}
with open(f"{out}/analysis/diversity_metrics.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        diversity[k] = v

# Top organism
top_organism = "unknown"
top_org_pct = "0"
with open(f"{out}/analysis/active_taxa.tsv") as f:
    next(f)
    first = f.readline().strip().split("\t")
    if len(first) >= 3:
        top_organism = first[0]
        top_org_pct = first[2]

# Top function
top_function = "unknown"
top_func_count = "0"
with open(f"{out}/analysis/functional_summary.tsv") as f:
    next(f)
    first = f.readline().strip().split("\t")
    if len(first) >= 2:
        top_function = first[0]
        top_func_count = first[1]

# Functional categories
func_cats = 0
with open(f"{out}/analysis/functional_enrichment.tsv") as f:
    next(f)
    for line in f:
        if line.strip():
            func_cats += 1

# Highly expressed genes
expressed_high = 0
with open(f"{out}/analysis/expression_cpm.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3 and float(parts[2]) > 100:
            expressed_high += 1

with open(f"{results}/report.csv", "w") as f:
    f.write("metric,value\n")
    f.write(f"total_reads_before,${READS_BEFORE}\n")
    f.write(f"total_reads_after,${READS_AFTER}\n")
    f.write(f"q30_rate_before,${Q30_BEFORE}\n")
    f.write(f"q30_rate_after,${Q30_AFTER}\n")
    f.write(f"rrna_reads,${RRNA_READS}\n")
    f.write(f"rrna_pct,${RRNA_PCT}\n")
    f.write(f"nonrrna_reads,${NONRRNA_READS}\n")
    f.write(f"total_contigs,${TOTAL_CONTIGS}\n")
    f.write(f"assembly_length,${ASSEMBLY_LENGTH}\n")
    f.write(f"largest_contig,${LARGEST_CONTIG}\n")
    f.write(f"assembly_gc_pct,${ASSEMBLY_GC}\n")
    f.write(f"assembly_n50,${N50}\n")
    f.write(f"predicted_genes,${PREDICTED_GENES}\n")
    f.write(f"annotated_genes,${ANNOTATED_GENES}\n")
    f.write(f"mapped_reads,${MAPPED_READS}\n")
    f.write(f"mapping_pct,${MAPPING_PCT}\n")
    f.write(f"expressed_genes,${EXPRESSED_GENES}\n")
    f.write(f"highly_expressed_genes,{expressed_high}\n")
    f.write(f"shannon_diversity,{diversity.get('shannon_diversity','')}\n")
    f.write(f"simpson_diversity,{diversity.get('simpson_diversity','')}\n")
    f.write(f"observed_richness,{diversity.get('observed_richness','')}\n")
    f.write(f"chao1_estimate,{diversity.get('chao1_estimate','')}\n")
    f.write(f"top_organism,{top_organism}\n")
    f.write(f"top_organism_pct,{top_org_pct}\n")
    f.write(f"top_function,{top_function}\n")
    f.write(f"functional_categories,{func_cats}\n")

print("Report written to results/report.csv")
PYEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RESULTS}/report.csv"
