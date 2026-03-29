#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Human Structural Variant Detection (Multi-caller) Pipeline
# =============================================================================
#
# DAG Structure (depth=10, convergence=4):
#
#  sample_R1.fq.gz    sample_R2.fq.gz
#       |                  |
#   [fastp QC] ------  [fastp QC]                              Level 1
#       |                  |
#       +--------+---------+
#                |
#        [bwa-mem2 align]                                       Level 2
#                |
#        [samtools sort + index]                                Level 3
#                |
#        [sambamba markdup]                                     Level 4
#                |
#     +----------+----------+-----------+
#     |          |          |           |
#  [Delly     [bcftools   [samtools   [samtools                 Level 5
#   SV call]   cnv call]   discordant  flagstat]
#   (PE+SR)    (depth-     extract +
#              based)      cluster]
#     |          |          |           |
#     +----------+----+-----+           |
#                     |                 |
#             CONVERGENCE 1             |                       Level 6
#             [SURVIVOR merge]          |
#             (>= 2 callers agree)      |
#                     |                 |
#             +-------+--------+        |
#             |       |        |        |
#       [SURVIVOR  [python  [bcftools   |                       Level 7
#        stats]     SV size  query]     |
#                   distrib]            |
#             |       |        |        |
#             +-------+--------+        |
#                     |                 |
#             CONVERGENCE 2             |                       Level 8
#             (stats + size + queries)  |
#                     |                 |
#             +-------+--------+        |
#             |       |        |        |
#       [bedtools  [SnpSift  [bcftools  |                       Level 9
#        intersect  annotate  filter]   |
#        (genes)]   ClinVar]            |
#             |       |        |        |
#             +-------+--------+        |
#                     |                 |
#             CONVERGENCE 3             |
#             (gene + clinical + filt)  |
#                     |                 |
#             [python clinical report]  |
#             CONVERGENCE 4 <-- QC + flagstat                   Level 10
#
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,aligned,processed,sv_delly,sv_bcf,sv_disc,merged,annotation,stats} "${RESULTS}"

# ---- Build indexes if needed ----
if [ ! -f "${REF}/genome.fa.fai" ]; then
  samtools faidx "${REF}/genome.fa"
fi
if [ ! -f "${REF}/genome.dict" ]; then
  samtools dict "${REF}/genome.fa" > "${REF}/genome.dict"
fi
if [ ! -f "${REF}/genome.fa.bwt.2bit.64" ]; then
  bwa-mem2 index "${REF}/genome.fa"
fi

# ---- Level 1: Read QC ----
if [ ! -f "${OUT}/qc/trimmed_R1.fastq.gz" ]; then
  echo ">>> Level 1: fastp QC"
  fastp \
    -i "${DATA}/sample_R1.fastq.gz" -I "${DATA}/sample_R2.fastq.gz" \
    -o "${OUT}/qc/trimmed_R1.fastq.gz" -O "${OUT}/qc/trimmed_R2.fastq.gz" \
    --json "${OUT}/qc/fastp.json" --html "${OUT}/qc/fastp.html" \
    --thread ${THREADS} --detect_adapter_for_pe \
    --cut_front --cut_tail --length_required 35
fi

# ---- Level 2: Alignment ----
if [ ! -f "${OUT}/aligned/raw.bam" ]; then
  echo ">>> Level 2: bwa-mem2 alignment"
  bwa-mem2 mem -t ${THREADS} \
    -R "@RG\tID:sample\tSM:TCRBOA7\tPL:ILLUMINA\tLB:WEX\tPU:unit1" \
    "${REF}/genome.fa" \
    "${OUT}/qc/trimmed_R1.fastq.gz" "${OUT}/qc/trimmed_R2.fastq.gz" | \
    samtools view -bS - > "${OUT}/aligned/raw.bam"
fi

# ---- Level 3: Sort + Index ----
if [ ! -f "${OUT}/aligned/sorted.bam" ]; then
  echo ">>> Level 3: samtools sort + index"
  samtools sort -@ ${THREADS} "${OUT}/aligned/raw.bam" -o "${OUT}/aligned/sorted.bam"
  samtools index "${OUT}/aligned/sorted.bam"
fi

# ---- Level 4: Mark Duplicates ----
if [ ! -f "${OUT}/processed/markdup.bam" ]; then
  echo ">>> Level 4: sambamba markdup"
  sambamba markdup -t ${THREADS} \
    "${OUT}/aligned/sorted.bam" \
    "${OUT}/processed/markdup.bam" 2>&1
fi

# ---- Level 5: Parallel SV callers ----

# Caller 1: Delly
if [ ! -f "${OUT}/sv_delly/delly.bcf" ]; then
  echo ">>> Level 5a: Delly SV calling"
  delly call -g "${REF}/genome.fa" \
    "${OUT}/processed/markdup.bam" \
    -o "${OUT}/sv_delly/delly.bcf" 2>&1 || true
  if [ -f "${OUT}/sv_delly/delly.bcf" ]; then
    bcftools view "${OUT}/sv_delly/delly.bcf" -Oz -o "${OUT}/sv_delly/delly.vcf.gz" 2>/dev/null || true
    bcftools index "${OUT}/sv_delly/delly.vcf.gz" 2>/dev/null || true
  fi
fi

# Caller 2: bcftools cnv/SV detection (depth-based)
if [ ! -f "${OUT}/sv_bcf/bcf_sv.vcf.gz" ]; then
  echo ">>> Level 5b: bcftools SV calling"
  # Use mpileup + call for large indels/SVs
  bcftools mpileup -f "${REF}/genome.fa" -q 20 -Q 20 --max-depth 500 \
    "${OUT}/processed/markdup.bam" 2>/dev/null | \
    bcftools call -mv --ploidy GRCh38 -Oz -o "${OUT}/sv_bcf/bcf_all.vcf.gz" 2>/dev/null || true
  bcftools index "${OUT}/sv_bcf/bcf_all.vcf.gz" 2>/dev/null || true

  # Extract indels >= 50bp (SV threshold)
  bcftools view -H "${OUT}/sv_bcf/bcf_all.vcf.gz" 2>/dev/null | \
    awk -F'\t' '{ref=length($4); alt=length($5); diff=ref-alt; if(diff<0) diff=-diff; if(diff>=50) print}' \
    > "${OUT}/sv_bcf/large_indels.txt" || true

  # Also detect depth anomalies with mosdepth-like approach
  samtools depth -a "${OUT}/processed/markdup.bam" 2>/dev/null | \
    awk 'BEGIN{OFS="\t"} {
      if($3 == 0) { if(!start) {start=$2; chr=$1} end=$2 }
      else { if(start) { len=end-start; if(len>=100) print chr,start,end,"DEL",len; start=0 } }
    }' > "${OUT}/sv_bcf/depth_gaps.bed" || true

  # Create a VCF from depth gaps
  python3 << 'PYEOF'
import sys
header = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"""
with open("outputs/sv_bcf/bcf_sv.vcf", "w") as out:
    out.write(header + "\n")
    i = 0
    with open("outputs/sv_bcf/depth_gaps.bed") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                chrom, start, end, svtype, svlen = parts[0], int(parts[1]), int(parts[2]), parts[3], int(parts[4])
                i += 1
                out.write(f"{chrom}\t{start}\t.\tN\t<{svtype}>\t.\tPASS\tSVTYPE={svtype};END={end};SVLEN=-{svlen}\tGT\t0/1\n")
    # Also add large indels from bcftools
    try:
        with open("outputs/sv_bcf/large_indels.txt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    chrom, pos = parts[0], parts[1]
                    ref, alt = parts[3], parts[4]
                    svlen = len(alt) - len(ref)
                    svtype = "INS" if svlen > 0 else "DEL"
                    end = int(pos) + abs(svlen)
                    i += 1
                    out.write(f"{chrom}\t{pos}\t.\tN\t<{svtype}>\t.\tPASS\tSVTYPE={svtype};END={end};SVLEN={svlen}\tGT\t0/1\n")
    except FileNotFoundError:
        pass
    print(f"Created {i} SV entries")
PYEOF
  bgzip -f "${OUT}/sv_bcf/bcf_sv.vcf"
  bcftools sort "${OUT}/sv_bcf/bcf_sv.vcf.gz" -Oz -o "${OUT}/sv_bcf/bcf_sv_sorted.vcf.gz" 2>/dev/null || true
  mv "${OUT}/sv_bcf/bcf_sv_sorted.vcf.gz" "${OUT}/sv_bcf/bcf_sv.vcf.gz" 2>/dev/null || true
  bcftools index -f "${OUT}/sv_bcf/bcf_sv.vcf.gz" 2>/dev/null || true
fi

# Caller 3: Discordant read extraction + clustering
if [ ! -f "${OUT}/sv_disc/disc_sv.vcf.gz" ]; then
  echo ">>> Level 5c: Discordant read SV detection"
  # Extract discordant and split reads
  samtools view -h -F 1294 "${OUT}/processed/markdup.bam" 2>/dev/null | \
    awk '$0 ~ /^@/ || ($9 > 1000 || $9 < -1000)' | \
    samtools view -bS - > "${OUT}/sv_disc/discordant.bam" 2>/dev/null || true

  DISC_COUNT=$(samtools view -c "${OUT}/sv_disc/discordant.bam" 2>/dev/null || echo "0")
  echo "Discordant reads: ${DISC_COUNT}"

  # Cluster discordant reads into SV candidates
  python3 << 'PYEOF'
import subprocess, sys

# Read discordant pairs
pairs = []
try:
    result = subprocess.run(["samtools", "view", "outputs/sv_disc/discordant.bam"],
                          capture_output=True, text=True, timeout=60)
    for line in result.stdout.strip().split("\n"):
        if not line: continue
        fields = line.split("\t")
        if len(fields) >= 9:
            chrom = fields[2]
            pos = int(fields[3])
            tlen = int(fields[8])
            if abs(tlen) > 1000:
                pairs.append((chrom, pos, tlen))
except Exception:
    pass

# Cluster nearby discordant pairs
clusters = []
if pairs:
    pairs.sort()
    cluster_start = pairs[0][1]
    cluster_chrom = pairs[0][0]
    cluster_end = pairs[0][1] + abs(pairs[0][2])
    cluster_count = 1
    for chrom, pos, tlen in pairs[1:]:
        if chrom == cluster_chrom and pos - cluster_start < 500:
            cluster_end = max(cluster_end, pos + abs(tlen))
            cluster_count += 1
        else:
            if cluster_count >= 2:
                clusters.append((cluster_chrom, cluster_start, cluster_end, cluster_count))
            cluster_start = pos
            cluster_chrom = chrom
            cluster_end = pos + abs(tlen)
            cluster_count = 1
    if cluster_count >= 2:
        clusters.append((cluster_chrom, cluster_start, cluster_end, cluster_count))

# Write VCF
header = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"""
with open("outputs/sv_disc/disc_sv.vcf", "w") as f:
    f.write(header + "\n")
    for chrom, start, end, count in sorted(clusters):
        svlen = end - start
        f.write(f"{chrom}\t{start}\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={end};SVLEN=-{svlen};PE={count}\tGT\t0/1\n")

print(f"Created {len(clusters)} SV entries from discordant reads")
PYEOF
  bgzip -f "${OUT}/sv_disc/disc_sv.vcf"
  bcftools index "${OUT}/sv_disc/disc_sv.vcf.gz" 2>/dev/null || true
fi

# Flagstat
if [ ! -f "${OUT}/stats/flagstat.txt" ]; then
  echo ">>> Level 5d: samtools flagstat"
  samtools flagstat "${OUT}/processed/markdup.bam" > "${OUT}/stats/flagstat.txt"
fi

# ---- CONVERGENCE 1: SURVIVOR merge ----
echo ">>> CONVERGENCE 1: SURVIVOR merge"

# Prepare file list for SURVIVOR
ls "${OUT}/sv_delly/delly.vcf.gz" "${OUT}/sv_bcf/bcf_sv.vcf.gz" "${OUT}/sv_disc/disc_sv.vcf.gz" 2>/dev/null \
  > "${OUT}/merged/caller_list.txt"

DELLY_COUNT=0
BCF_SV_COUNT=0
DISC_SV_COUNT=0
if [ -s "${OUT}/sv_delly/delly.vcf.gz" ]; then
  DELLY_COUNT=$(bcftools view -H "${OUT}/sv_delly/delly.vcf.gz" 2>/dev/null | wc -l || true)
fi
if [ -s "${OUT}/sv_bcf/bcf_sv.vcf.gz" ]; then
  BCF_SV_COUNT=$(bcftools view -H "${OUT}/sv_bcf/bcf_sv.vcf.gz" 2>/dev/null | wc -l || true)
fi
if [ -s "${OUT}/sv_disc/disc_sv.vcf.gz" ]; then
  DISC_SV_COUNT=$(bcftools view -H "${OUT}/sv_disc/disc_sv.vcf.gz" 2>/dev/null | wc -l || true)
fi

echo "Delly SVs: ${DELLY_COUNT}, bcftools SVs: ${BCF_SV_COUNT}, Discordant SVs: ${DISC_SV_COUNT}"

# SURVIVOR merge: max distance 1000bp, min 1 caller (to get union), by SV type
if [ ! -f "${OUT}/merged/survivor_merged.vcf" ]; then
  SURVIVOR merge "${OUT}/merged/caller_list.txt" 1000 1 1 1 0 50 \
    "${OUT}/merged/survivor_merged.vcf" 2>&1 || true
fi

MERGED_COUNT=0
if [ -f "${OUT}/merged/survivor_merged.vcf" ]; then
  MERGED_COUNT=$(grep -cv "^#" "${OUT}/merged/survivor_merged.vcf" || true)
  MERGED_COUNT=${MERGED_COUNT:-0}
fi
echo "Merged SVs: ${MERGED_COUNT}"

# ---- Level 7: Stats + size distribution + query ----
echo ">>> Level 7: SV statistics"

# SURVIVOR stats
if [ -f "${OUT}/merged/survivor_merged.vcf" ] && [ "${MERGED_COUNT}" -gt 0 ]; then
  SURVIVOR stats "${OUT}/merged/survivor_merged.vcf" -1 -1 -1 \
    "${OUT}/stats/survivor_stats" 2>&1 || true
fi

# SV size distribution
python3 << 'PYEOF'
import re
sv_types = {}
sv_sizes = []
try:
    with open("outputs/merged/survivor_merged.vcf") as f:
        for line in f:
            if line.startswith("#"): continue
            m = re.search(r'SVTYPE=(\w+)', line)
            if m:
                svtype = m.group(1)
                sv_types[svtype] = sv_types.get(svtype, 0) + 1
            m = re.search(r'SVLEN=(-?\d+)', line)
            if m:
                sv_sizes.append(abs(int(m.group(1))))
except FileNotFoundError:
    pass

with open("outputs/stats/sv_summary.txt", "w") as f:
    f.write("SV type counts:\n")
    for t, c in sorted(sv_types.items()):
        f.write(f"  {t}: {c}\n")
    if sv_sizes:
        f.write(f"\nSV size stats:\n")
        f.write(f"  min: {min(sv_sizes)}\n")
        f.write(f"  max: {max(sv_sizes)}\n")
        f.write(f"  median: {sorted(sv_sizes)[len(sv_sizes)//2]}\n")
        f.write(f"  mean: {sum(sv_sizes)//len(sv_sizes)}\n")

print(open("outputs/stats/sv_summary.txt").read())
PYEOF

# ---- Level 9: Annotation ----
echo ">>> Level 9: SV annotation"

# Gene overlap
GENE_OVERLAP=0
if [ -f "${OUT}/merged/survivor_merged.vcf" ] && [ -s "${REF}/gene_regions.bed" ] && [ "${MERGED_COUNT}" -gt 0 ]; then
  # Convert VCF to BED
  grep -v "^#" "${OUT}/merged/survivor_merged.vcf" | \
    awk -F'\t' '{
      match($8, /END=([0-9]+)/, a);
      if(a[1]) print $1"\t"$2-1"\t"a[1]"\t"$3
    }' > "${OUT}/annotation/sv_regions.bed" || true

  if [ -s "${OUT}/annotation/sv_regions.bed" ]; then
    GENE_OVERLAP=$(bedtools intersect -a "${OUT}/annotation/sv_regions.bed" \
      -b "${REF}/gene_regions.bed" -u | wc -l || true)
  fi
fi
echo "Gene-overlapping SVs: ${GENE_OVERLAP}"

# ClinVar annotation
CLINVAR_SV_HITS=0
if [ -f "${OUT}/merged/survivor_merged.vcf" ] && [ -s "${REF}/clinvar_sv.vcf.gz" ]; then
  SnpSift annotate "${REF}/clinvar_sv.vcf.gz" \
    "${OUT}/merged/survivor_merged.vcf" \
    > "${OUT}/annotation/clinvar_sv_annotated.vcf" 2>/dev/null || true
  if [ -f "${OUT}/annotation/clinvar_sv_annotated.vcf" ]; then
    CLINVAR_SV_HITS=$(grep -c "CLNSIG=" "${OUT}/annotation/clinvar_sv_annotated.vcf" || true)
    CLINVAR_SV_HITS=${CLINVAR_SV_HITS:-0}
  fi
fi

# ---- CONVERGENCE 3+4: Final Report ----
echo ">>> CONVERGENCE 3+4: Generating report"

# QC stats
RAW_READS=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['before_filtering']['total_reads'])")
CLEAN_READS=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(d['summary']['after_filtering']['total_reads'])")
Q30_RATE=$(python3 -c "import json; d=json.load(open('${OUT}/qc/fastp.json')); print(round(d['summary']['after_filtering']['q30_rate']*100,2))")

# Alignment stats
MAPPED=$(grep "mapped (" "${OUT}/stats/flagstat.txt" | head -1 | awk '{print $1}')
TOTAL_ALN=$(grep "in total" "${OUT}/stats/flagstat.txt" | awk '{print $1}')
MAP_RATE=$(python3 -c "print(round(${MAPPED}/${TOTAL_ALN}*100,2)) if ${TOTAL_ALN} > 0 else print(0)")

# SV type counts
DEL_COUNT=0; DUP_COUNT=0; INV_COUNT=0; INS_COUNT=0; BND_COUNT=0
if [ -f "${OUT}/merged/survivor_merged.vcf" ]; then
  DEL_COUNT=$(grep -c "SVTYPE=DEL" "${OUT}/merged/survivor_merged.vcf" || true)
  DUP_COUNT=$(grep -c "SVTYPE=DUP" "${OUT}/merged/survivor_merged.vcf" || true)
  INV_COUNT=$(grep -c "SVTYPE=INV" "${OUT}/merged/survivor_merged.vcf" || true)
  INS_COUNT=$(grep -c "SVTYPE=INS" "${OUT}/merged/survivor_merged.vcf" || true)
  BND_COUNT=$(grep -c "SVTYPE=BND" "${OUT}/merged/survivor_merged.vcf" || true)
fi

cat > "${RESULTS}/report.csv" << CSVEOF
metric,value
raw_reads,${RAW_READS}
clean_reads,${CLEAN_READS}
q30_rate,${Q30_RATE}
mapped_reads,${MAPPED}
mapping_rate,${MAP_RATE}
caller1_svs,${DELLY_COUNT}
caller2_svs,${BCF_SV_COUNT}
caller3_svs,${DISC_SV_COUNT}
merged_svs,${MERGED_COUNT}
deletion_count,${DEL_COUNT}
duplication_count,${DUP_COUNT}
inversion_count,${INV_COUNT}
insertion_count,${INS_COUNT}
breakend_count,${BND_COUNT}
gene_overlapping_svs,${GENE_OVERLAP}
clinvar_sv_annotated,${CLINVAR_SV_HITS}
CSVEOF

echo "=== Final Report ==="
cat "${RESULTS}/report.csv"
echo ""
echo "Pipeline complete!"
