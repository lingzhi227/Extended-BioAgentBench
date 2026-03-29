#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Pharmacogenomics Pipeline: CYP2D6 Star Allele Calling
# ============================================================
# DAG structure (depth=9, convergence=4):
#
#  sample.R1.fq.gz   sample.R2.fq.gz
#        │                │
#    [fastp QC]        [fastp QC]                     Level 1
#        │                │
#        └────────┬───────┘
#                 │
#         [bwa-mem2 align to chr22]                   Level 2
#                 │
#         [samtools sort + index]                     Level 3
#                 │
#         [picard MarkDuplicates]                     Level 4
#                 │
#         ┌───────┼───────────────────┐
#         │       │                   │
#    [gatk HC    [samtools           [mosdepth         Level 5
#     CYP region  view CYP2D6]       CYP2D6
#     intervals]       │              coverage]
#         │       [samtools depth]        │
#         │            │                  │
#    [bcftools  [CONVERGENCE 1] ◄─────────┘           Level 6
#     norm]     (CYP2D6 BAM + coverage)
#         │            │
#         │     [python CN estimation]                Level 7
#         │     (depth ratio CYP2D6/CYP2D7)
#         │            │
#         └────────────┤
#                      │
#              [CONVERGENCE 2]                        Level 8
#              (VCF + CN + coverage)
#                      │
#         ┌────────────┼────────────┐
#         │            │            │
#   [python star  [bcftools     [python               Level 8
#    allele        annotate      drug-gene
#    calling]      (rsID)]       interaction
#         │            │          lookup]
#         │            │            │
#         └────────────┼────────────┘
#                      │
#              [CONVERGENCE 3]                        Level 8
#              (diplotype + annotation + drugs)
#                      │
#         ┌────────────┼────────┐
#         │            │        │
#   [python       [python    [python                  Level 9
#    metabolizer   activity   phenotype
#    status]       score]     prediction]
#         │            │        │
#         └────────────┼────────┘
#                      │
#              [CONVERGENCE 4]
#              [clinical report]
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,align,dedup,variants,coverage,cyp2d6,star_alleles,drugs,report}
mkdir -p "${RESULTS}"

R1="${DATA}/HG002_chr22_R1.fastq.gz"
R2="${DATA}/HG002_chr22_R2.fastq.gz"
REFFA="${REF}/chr22.fasta"

# CYP2D6 region (GRCh38)
CYP2D6_REGION="chr22:42126000-42132500"
CYP2D7_REGION="chr22:42137000-42143500"
BROAD_PGX_REGION="chr22:42120000-42150000"

# ============================================================
# Level 1: Read QC with fastp
# ============================================================
echo "=== Level 1: fastp QC ==="
if [ ! -f "${OUT}/qc/fastp.json" ]; then
  fastp -i "${R1}" -I "${R2}" \
    -o "${OUT}/qc/trimmed_R1.fastq.gz" \
    -O "${OUT}/qc/trimmed_R2.fastq.gz" \
    -j "${OUT}/qc/fastp.json" \
    -h "${OUT}/qc/fastp.html" \
    -w ${THREADS} \
    > "${OUT}/qc/fastp.log" 2>&1
  echo "  fastp done"
fi

# ============================================================
# Level 2: Alignment with bwa-mem2
# ============================================================
echo "=== Level 2: bwa-mem2 alignment ==="
if [ ! -f "${OUT}/align/aligned.bam" ]; then
  bwa-mem2 mem -t ${THREADS} \
    -R "@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tLB:lib1" \
    "${REFFA}" \
    "${OUT}/qc/trimmed_R1.fastq.gz" \
    "${OUT}/qc/trimmed_R2.fastq.gz" 2>/dev/null | \
    samtools view -bS -@ ${THREADS} - > "${OUT}/align/aligned.bam" 2>/dev/null
  echo "  Alignment done"
fi

# ============================================================
# Level 3: Sort + index
# ============================================================
echo "=== Level 3: Sort + index ==="
if [ ! -f "${OUT}/align/sorted.bam" ]; then
  samtools sort -@ ${THREADS} "${OUT}/align/aligned.bam" -o "${OUT}/align/sorted.bam"
  samtools index "${OUT}/align/sorted.bam"
  echo "  Sort done"
fi

# ============================================================
# Level 4: Mark duplicates
# ============================================================
echo "=== Level 4: Mark duplicates ==="
if [ ! -f "${OUT}/dedup/dedup.bam" ]; then
  picard MarkDuplicates \
    I="${OUT}/align/sorted.bam" \
    O="${OUT}/dedup/dedup.bam" \
    M="${OUT}/dedup/metrics.txt" \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT \
    > "${OUT}/dedup/picard.log" 2>&1
  samtools index "${OUT}/dedup/dedup.bam"
  echo "  Dedup done"
fi

# ============================================================
# Level 5: Three parallel branches from dedup BAM
# ============================================================
echo "=== Level 5: Variant calling + CYP2D6 extraction + coverage ==="

# Branch 5a: GATK HaplotypeCaller on PGx region
if [ ! -f "${OUT}/variants/raw.vcf.gz" ]; then
  # Create intervals file
  echo "${BROAD_PGX_REGION}" | tr ':' '\t' | tr '-' '\t' > "${OUT}/variants/pgx_region.bed"

  gatk HaplotypeCaller \
    -R "${REFFA}" \
    -I "${OUT}/dedup/dedup.bam" \
    -L "${BROAD_PGX_REGION}" \
    -O "${OUT}/variants/raw.vcf.gz" \
    --native-pair-hmm-threads ${THREADS} \
    > "${OUT}/variants/gatk.log" 2>&1
  echo "  GATK HC done"
fi

# Normalize VCF
if [ ! -f "${OUT}/variants/norm.vcf.gz" ]; then
  bcftools norm -f "${REFFA}" \
    -m -both \
    "${OUT}/variants/raw.vcf.gz" | \
    bcftools view -i 'QUAL>=30' -Oz -o "${OUT}/variants/norm.vcf.gz"
  bcftools index "${OUT}/variants/norm.vcf.gz"
  echo "  bcftools norm done"
fi

# Branch 5b: Extract CYP2D6 region reads
if [ ! -f "${OUT}/cyp2d6/cyp2d6_reads.bam" ]; then
  samtools view -b "${OUT}/dedup/dedup.bam" ${CYP2D6_REGION} > "${OUT}/cyp2d6/cyp2d6_reads.bam"
  samtools index "${OUT}/cyp2d6/cyp2d6_reads.bam"

  # Also extract CYP2D7 region for CN estimation
  samtools view -b "${OUT}/dedup/dedup.bam" ${CYP2D7_REGION} > "${OUT}/cyp2d6/cyp2d7_reads.bam"
  samtools index "${OUT}/cyp2d6/cyp2d7_reads.bam"
  echo "  CYP2D6/D7 extraction done"
fi

# Per-base depth
if [ ! -f "${OUT}/cyp2d6/cyp2d6_depth.txt" ]; then
  samtools depth -r ${CYP2D6_REGION} "${OUT}/dedup/dedup.bam" > "${OUT}/cyp2d6/cyp2d6_depth.txt"
  samtools depth -r ${CYP2D7_REGION} "${OUT}/dedup/dedup.bam" > "${OUT}/cyp2d6/cyp2d7_depth.txt"
  echo "  Per-base depth done"
fi

# Branch 5c: mosdepth coverage
if [ ! -f "${OUT}/coverage/pgx_region.mosdepth.summary.txt" ]; then
  echo -e "chr22\t42126000\t42132500\tCYP2D6" > "${OUT}/coverage/pgx_regions.bed"
  echo -e "chr22\t42137000\t42143500\tCYP2D7" >> "${OUT}/coverage/pgx_regions.bed"

  mosdepth --by "${OUT}/coverage/pgx_regions.bed" \
    --no-per-base \
    -t ${THREADS} \
    "${OUT}/coverage/pgx_region" \
    "${OUT}/dedup/dedup.bam" \
    > /dev/null 2>&1
  echo "  mosdepth done"
fi

# ============================================================
# Level 6: CONVERGENCE 1 — CYP2D6 BAM + coverage
# ============================================================
echo "=== Level 6: Convergence 1 ==="

# ============================================================
# Level 7: Copy number estimation
# ============================================================
echo "=== Level 7: CN estimation ==="
if [ ! -f "${OUT}/cyp2d6/cn_estimate.json" ]; then
  python3 << 'PYEOF'
import json
import statistics

# Read CYP2D6 depth
cyp2d6_depths = []
with open("outputs/cyp2d6/cyp2d6_depth.txt") as f:
    for line in f:
        parts = line.strip().split('\t')
        cyp2d6_depths.append(int(parts[2]))

# Read CYP2D7 depth
cyp2d7_depths = []
with open("outputs/cyp2d6/cyp2d7_depth.txt") as f:
    for line in f:
        parts = line.strip().split('\t')
        cyp2d7_depths.append(int(parts[2]))

cyp2d6_mean = statistics.mean(cyp2d6_depths) if cyp2d6_depths else 0
cyp2d7_mean = statistics.mean(cyp2d7_depths) if cyp2d7_depths else 0

# CN estimation: CYP2D7 is always diploid (2 copies), use as reference
if cyp2d7_mean > 0:
    cn_ratio = cyp2d6_mean / cyp2d7_mean
    cn_estimate = round(cn_ratio * 2)  # Diploid reference = 2
else:
    cn_ratio = 1.0
    cn_estimate = 2

result = {
    "cyp2d6_mean_depth": round(cyp2d6_mean, 1),
    "cyp2d7_mean_depth": round(cyp2d7_mean, 1),
    "depth_ratio": round(cn_ratio, 3),
    "copy_number_estimate": cn_estimate,
    "has_deletion": cn_estimate < 2,
    "has_duplication": cn_estimate > 2
}

with open("outputs/cyp2d6/cn_estimate.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"CYP2D6 mean depth: {cyp2d6_mean:.1f}")
print(f"CYP2D7 mean depth: {cyp2d7_mean:.1f}")
print(f"Ratio: {cn_ratio:.3f}, CN estimate: {cn_estimate}")
PYEOF
fi

# ============================================================
# Level 8: CONVERGENCE 2 — VCF + CN + coverage
# ============================================================
echo "=== Level 8: Convergence 2 + Star allele calling ==="

# Branch 8a: Star allele calling
if [ ! -f "${OUT}/star_alleles/diplotype.json" ]; then
  python3 << 'PYEOF'
import json
import subprocess

# Load PharmCAT CYP2D6 definitions
with open("reference/CYP2D6_translation.json") as f:
    pharmcat = json.load(f)

# Parse VCF for CYP2D6 region variants
variants = {}
result = subprocess.run(
    ["bcftools", "query", "-f", "%POS\t%REF\t%ALT\t[%GT]\n",
     "-r", "chr22:42126000-42132500",
     "outputs/variants/norm.vcf.gz"],
    capture_output=True, text=True
)

for line in result.stdout.strip().split('\n'):
    if not line:
        continue
    parts = line.split('\t')
    pos = int(parts[0])
    ref = parts[1]
    alt = parts[2]
    gt = parts[3]
    variants[pos] = {"ref": ref, "alt": alt, "gt": gt}

# Get the variant positions from PharmCAT definition
star_variants = {}
if "variants" in pharmcat:
    for v in pharmcat["variants"]:
        pos = v.get("position")
        if pos:
            star_variants[int(pos)] = v

# Simple star allele matching
# Check which defining variants are present
defining_snps = []
for pos, info in sorted(variants.items()):
    gt = info["gt"]
    if gt in ["0/1", "1/1", "0|1", "1|1"]:
        defining_snps.append({
            "position": pos,
            "ref": info["ref"],
            "alt": info["alt"],
            "genotype": gt,
            "zygosity": "heterozygous" if "0" in gt else "homozygous"
        })

# Load CN estimate
with open("outputs/cyp2d6/cn_estimate.json") as f:
    cn = json.load(f)

# Basic diplotype assignment based on CN and variants
# *1 = reference (normal function), *5 = deletion, *xN = duplication
if cn["has_deletion"]:
    allele1 = "*5"  # deletion
    allele2 = "*1"
elif cn["has_duplication"]:
    allele1 = "*1"
    allele2 = "*1xN"
else:
    # Check known defining variants
    # Most common star alleles by defining SNP:
    # *2: rs16947 (42130692 G>A) — normal function
    # *4: rs3892097 (42128945 G>A) — no function
    # *10: rs1065852 (42130728 C>T) — decreased function
    # *41: rs28371725 (42127803 C>T) — decreased function
    allele1 = "*1"
    allele2 = "*1"

    known_stars = {
        42130692: ("*2", "normal_function"),
        42128945: ("*4", "no_function"),
        42130728: ("*10", "decreased_function"),
        42127803: ("*41", "decreased_function"),
        42126611: ("*3", "no_function"),
        42129132: ("*6", "no_function"),
        42129770: ("*9", "decreased_function"),
        42129819: ("*17", "decreased_function"),
    }

    found_alleles = []
    for snp in defining_snps:
        pos = snp["position"]
        if pos in known_stars:
            star, function = known_stars[pos]
            found_alleles.append({"star": star, "function": function, "zygosity": snp["zygosity"]})

    if found_alleles:
        if found_alleles[0]["zygosity"] == "homozygous":
            allele1 = found_alleles[0]["star"]
            allele2 = found_alleles[0]["star"]
        else:
            allele1 = found_alleles[0]["star"]
            allele2 = found_alleles[1]["star"] if len(found_alleles) > 1 else "*1"

diplotype = f"{allele1}/{allele2}"

result = {
    "diplotype": diplotype,
    "allele1": allele1,
    "allele2": allele2,
    "copy_number": cn["copy_number_estimate"],
    "defining_variants_found": len(defining_snps),
    "total_variants_in_region": len(variants),
    "known_star_allele_variants": [
        {"star": a["star"], "function": a["function"]}
        for a in (found_alleles if 'found_alleles' in dir() else [])
    ]
}

with open("outputs/star_alleles/diplotype.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"Diplotype: {diplotype}")
print(f"Defining variants: {len(defining_snps)}")
PYEOF
fi

# Branch 8b: Annotate VCF with rsIDs
if [ ! -f "${OUT}/star_alleles/annotated_variants.tsv" ]; then
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\n' \
    -r chr22:42126000-42132500 \
    "${OUT}/variants/norm.vcf.gz" > "${OUT}/star_alleles/annotated_variants.tsv" 2>/dev/null
  echo "  Variant annotation done"
fi

# Branch 8c: Drug-gene interaction lookup
if [ ! -f "${OUT}/drugs/drug_interactions.json" ]; then
  python3 << 'PYEOF'
import json

# Load diplotype
with open("outputs/star_alleles/diplotype.json") as f:
    dp = json.load(f)

# Load CPIC pairs
with open("reference/cpic_cyp2d6_pairs.json") as f:
    pairs = json.load(f)

# Load allele function data
with open("reference/cpic_cyp2d6_alleles.json") as f:
    alleles = json.load(f)

# Build function lookup
allele_functions = {}
for a in alleles:
    name = a.get("name", "")
    func = a.get("functionalstatus", a.get("clinicalfunctionalstatus", ""))
    if name and func:
        allele_functions[name] = func

# Activity score calculation
activity_scores = {
    "normal_function": 1.0,
    "Normal function": 1.0,
    "decreased_function": 0.5,
    "Decreased function": 0.5,
    "no_function": 0.0,
    "No function": 0.0,
    "uncertain_function": 0.5,
    "Uncertain function": 0.5,
}

a1_func = allele_functions.get(dp["allele1"], "normal_function")
a2_func = allele_functions.get(dp["allele2"], "normal_function")
a1_score = activity_scores.get(a1_func, 1.0)
a2_score = activity_scores.get(a2_func, 1.0)
total_score = a1_score + a2_score

# Phenotype determination based on activity score
if total_score >= 2.25:
    phenotype = "Ultrarapid Metabolizer"
elif total_score >= 1.25:
    phenotype = "Normal Metabolizer"
elif total_score >= 0.25:
    phenotype = "Intermediate Metabolizer"
else:
    phenotype = "Poor Metabolizer"

# Key affected drugs
key_drugs = []
for p in pairs:
    drug = p.get("drugname", p.get("drugid", "Unknown"))
    level = p.get("cpiclevel", "")
    if level in ["A", "B"]:
        key_drugs.append({"drug": drug, "cpic_level": level})

result = {
    "diplotype": dp["diplotype"],
    "allele1_function": a1_func,
    "allele2_function": a2_func,
    "activity_score": total_score,
    "phenotype": phenotype,
    "total_cpic_drug_pairs": len(pairs),
    "level_a_b_drugs": len(key_drugs),
    "key_affected_drugs": key_drugs[:10]  # top 10
}

with open("outputs/drugs/drug_interactions.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"Phenotype: {phenotype} (activity score: {total_score})")
print(f"Key drugs: {len(key_drugs)} CPIC Level A/B")
PYEOF
fi

# ============================================================
# CONVERGENCE 3 + Level 9: Clinical report
# ============================================================
echo "=== Level 9: Clinical report ==="

python3 << 'PYEOF'
import json
import os

os.chdir(os.path.dirname(os.path.abspath("run_script.sh")) or ".")

metrics = {}

# QC metrics
with open("outputs/qc/fastp.json") as f:
    qc = json.load(f)
    metrics["total_reads_before_qc"] = qc["summary"]["before_filtering"]["total_reads"]
    metrics["total_reads_after_qc"] = qc["summary"]["after_filtering"]["total_reads"]
    metrics["q30_rate_pct"] = round(qc["summary"]["after_filtering"]["q30_rate"] * 100, 2)

# Alignment stats
flagstat_lines = os.popen("samtools flagstat outputs/dedup/dedup.bam 2>/dev/null").readlines()
for line in flagstat_lines:
    if "mapped (" in line:
        parts = line.strip().split()
        metrics["mapped_reads"] = int(parts[0])
        import re
        pct = re.search(r'\(([\d.]+)%', line)
        if pct:
            metrics["mapping_rate_pct"] = float(pct.group(1))
        break

# Dedup metrics
with open("outputs/dedup/metrics.txt") as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("## METRICS"):
            header = lines[i+1].strip().split('\t')
            values = lines[i+2].strip().split('\t')
            d = dict(zip(header, values))
            metrics["duplicate_rate_pct"] = round(float(d.get("PERCENT_DUPLICATION", 0)) * 100, 2)
            break

# CYP2D6 coverage
with open("outputs/cyp2d6/cn_estimate.json") as f:
    cn = json.load(f)
    metrics["cyp2d6_mean_depth"] = cn["cyp2d6_mean_depth"]
    metrics["cyp2d7_mean_depth"] = cn["cyp2d7_mean_depth"]
    metrics["cyp2d6_cyp2d7_depth_ratio"] = cn["depth_ratio"]
    metrics["copy_number_estimate"] = cn["copy_number_estimate"]

# Variant stats
import subprocess
result = subprocess.run(
    ["bcftools", "stats", "-r", "chr22:42126000-42132500",
     "outputs/variants/norm.vcf.gz"],
    capture_output=True, text=True
)
for line in result.stdout.split('\n'):
    if line.startswith("SN") and "number of records:" in line:
        metrics["variants_in_cyp2d6_region"] = int(line.strip().split('\t')[-1])
    if line.startswith("SN") and "number of SNPs:" in line:
        metrics["snps_in_cyp2d6_region"] = int(line.strip().split('\t')[-1])

# Star allele results
with open("outputs/star_alleles/diplotype.json") as f:
    dp = json.load(f)
    metrics["diplotype"] = dp["diplotype"]
    metrics["defining_variants_found"] = dp["defining_variants_found"]

# Drug interaction results
with open("outputs/drugs/drug_interactions.json") as f:
    drugs = json.load(f)
    metrics["activity_score"] = drugs["activity_score"]
    metrics["metabolizer_phenotype"] = drugs["phenotype"]
    metrics["cpic_drug_pairs_total"] = drugs["total_cpic_drug_pairs"]
    metrics["cpic_level_ab_drugs"] = drugs["level_a_b_drugs"]

# mosdepth summary
mosdepth_file = "outputs/coverage/pgx_region.regions.bed.gz"
if os.path.exists(mosdepth_file):
    import gzip
    with gzip.open(mosdepth_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                region_name = parts[3] if len(parts) > 3 else parts[0]
                mean_cov = float(parts[4]) if len(parts) > 4 else float(parts[3])
                if "CYP2D6" in line and "CYP2D7" not in line:
                    metrics["mosdepth_cyp2d6_coverage"] = round(mean_cov, 1)
                elif "CYP2D7" in line:
                    metrics["mosdepth_cyp2d7_coverage"] = round(mean_cov, 1)

# Write report
with open("results/report.csv", 'w') as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("=== Report ===")
for k, v in metrics.items():
    print(f"  {k} = {v}")
PYEOF

echo "=== Pipeline complete ==="
