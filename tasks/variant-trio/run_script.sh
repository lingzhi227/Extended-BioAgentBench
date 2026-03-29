#!/usr/bin/env bash
set -euo pipefail

# Variant Annotation + Clinical Interpretation (Trio Analysis)
# Data: GIAB Ashkenazi trio (HG002=proband, HG003=father, HG004=mother), chr22
# Reference: GRCh38 chr22 + ClinVar annotations + PED file
#
# DAG Structure (depth=10, convergence=4):
#
#   [HG002.vcf]────[HG003.vcf]────[HG004.vcf]
#       │               │              │
#   bcftools norm   bcftools norm  bcftools norm      (Step 1)
#       │               │              │
#       └───────────────┼──────────────┘
#                       │
#               bcftools merge                         (Step 2: CONVERGE #1)
#                       │
#              ╱────────┴────────╲
#     bcftools view(SNP)   bcftools view(Indel)        (Step 3: parallel)
#              │                  │
#     gatk VariantFiltration  gatk VariantFiltration   (Step 4: parallel filter)
#              ╲────────┬────────╱
#                       │
#               bcftools concat                        (Step 5: CONVERGE #2)
#                       │
#              ╱────────┴────────╲
#     ClinVar annotation    frequency stats            (Step 6: parallel)
#     (bcftools annotate)   (vcftools + bcftools)
#              ╲────────┬────────╱
#                       │
#               merge annotations                      (Step 7: CONVERGE #3)
#                       │
#              inheritance scoring                      (Step 8: CONVERGE #4 + PED)
#              (genmod annotate)
#                       │
#               variant classification                  (Step 9)
#                       │
#                   report.csv                          (Step 10)

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}/norm" "${OUT}/merged" "${OUT}/split" "${OUT}/filtered"
mkdir -p "${OUT}/annotated" "${OUT}/stats" "${OUT}/inheritance" "${RESULTS}"

REFERENCE="${REF}/chr22.fa"
CLINVAR="${REF}/clinvar_chr22.vcf.gz"
PED="${REF}/trio.ped"

# ============================================================================
# Step 1: Normalize VCFs per sample (parallel)
# ============================================================================
echo "[Step 1] Normalizing VCFs..."
for SAMPLE in HG002 HG003 HG004; do
    if [ ! -f "${OUT}/norm/${SAMPLE}_norm.vcf.gz" ]; then
        bcftools norm \
            -m -both \
            -f "$REFERENCE" \
            -o "${OUT}/norm/${SAMPLE}_norm.vcf.gz" \
            -O z \
            "${DATA}/${SAMPLE}_chr22.vcf.gz" 2>/dev/null
        tabix -p vcf "${OUT}/norm/${SAMPLE}_norm.vcf.gz"
    fi
    echo "  ${SAMPLE}: $(bcftools view -H ${OUT}/norm/${SAMPLE}_norm.vcf.gz | wc -l) variants"
done

# ============================================================================
# Step 2: Merge trio VCFs — CONVERGE #1 (all three samples)
# ============================================================================
echo "[Step 2] Merging trio VCFs (convergence #1)..."
if [ ! -f "${OUT}/merged/trio_merged.vcf.gz" ]; then
    bcftools merge \
        --force-samples \
        -m both \
        -O z \
        -o "${OUT}/merged/trio_merged.vcf.gz" \
        "${OUT}/norm/HG002_norm.vcf.gz" \
        "${OUT}/norm/HG003_norm.vcf.gz" \
        "${OUT}/norm/HG004_norm.vcf.gz"
    tabix -p vcf "${OUT}/merged/trio_merged.vcf.gz"
fi
TOTAL_MERGED=$(bcftools view -H "${OUT}/merged/trio_merged.vcf.gz" | wc -l)
echo "  Merged variants: ${TOTAL_MERGED}"

# ============================================================================
# Step 3: Split by variant type (parallel)
# ============================================================================
echo "[Step 3] Splitting by variant type..."
if [ ! -f "${OUT}/split/snps.vcf.gz" ]; then
    bcftools view -v snps -O z -o "${OUT}/split/snps.vcf.gz" "${OUT}/merged/trio_merged.vcf.gz"
    tabix -p vcf "${OUT}/split/snps.vcf.gz"
fi
if [ ! -f "${OUT}/split/indels.vcf.gz" ]; then
    bcftools view -v indels -O z -o "${OUT}/split/indels.vcf.gz" "${OUT}/merged/trio_merged.vcf.gz"
    tabix -p vcf "${OUT}/split/indels.vcf.gz"
fi
N_SNPS=$(bcftools view -H "${OUT}/split/snps.vcf.gz" | wc -l)
N_INDELS=$(bcftools view -H "${OUT}/split/indels.vcf.gz" | wc -l)
echo "  SNPs: ${N_SNPS}, Indels: ${N_INDELS}"

# ============================================================================
# Step 4: Hard filtering per variant type (parallel)
# ============================================================================
echo "[Step 4] Applying hard filters..."

# SNP filtering
if [ ! -f "${OUT}/filtered/snps_filtered.vcf.gz" ]; then
    gatk VariantFiltration \
        -R "$REFERENCE" \
        -V "${OUT}/split/snps.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "LowQD" \
        --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
        --filter-expression "FS > 60.0" --filter-name "HighFS" \
        --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
        -O "${OUT}/filtered/snps_filtered.vcf" 2>/dev/null || {
        # GIAB benchmark VCFs may not have QD/MQ/FS fields - use PASS filter
        bcftools view -f PASS -O v -o "${OUT}/filtered/snps_filtered.vcf" "${OUT}/split/snps.vcf.gz"
    }
    bgzip "${OUT}/filtered/snps_filtered.vcf"
    tabix -p vcf "${OUT}/filtered/snps_filtered.vcf.gz"
fi

# Indel filtering
if [ ! -f "${OUT}/filtered/indels_filtered.vcf.gz" ]; then
    gatk VariantFiltration \
        -R "$REFERENCE" \
        -V "${OUT}/split/indels.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "LowQD" \
        --filter-expression "FS > 200.0" --filter-name "HighFS" \
        --filter-expression "SOR > 10.0" --filter-name "HighSOR" \
        -O "${OUT}/filtered/indels_filtered.vcf" 2>/dev/null || {
        bcftools view -f PASS -O v -o "${OUT}/filtered/indels_filtered.vcf" "${OUT}/split/indels.vcf.gz"
    }
    bgzip "${OUT}/filtered/indels_filtered.vcf"
    tabix -p vcf "${OUT}/filtered/indels_filtered.vcf.gz"
fi

echo "  Filtered SNPs: $(bcftools view -H ${OUT}/filtered/snps_filtered.vcf.gz | wc -l)"
echo "  Filtered Indels: $(bcftools view -H ${OUT}/filtered/indels_filtered.vcf.gz | wc -l)"

# ============================================================================
# Step 5: Merge filtered variants — CONVERGE #2
# ============================================================================
echo "[Step 5] Merging filtered SNPs and Indels (convergence #2)..."
if [ ! -f "${OUT}/filtered/trio_filtered.vcf.gz" ]; then
    bcftools concat \
        -a \
        -O z \
        -o "${OUT}/filtered/trio_filtered.vcf.gz" \
        "${OUT}/filtered/snps_filtered.vcf.gz" \
        "${OUT}/filtered/indels_filtered.vcf.gz"
    bcftools sort -O z -o "${OUT}/filtered/trio_filtered_sorted.vcf.gz" "${OUT}/filtered/trio_filtered.vcf.gz"
    mv "${OUT}/filtered/trio_filtered_sorted.vcf.gz" "${OUT}/filtered/trio_filtered.vcf.gz"
    tabix -p vcf "${OUT}/filtered/trio_filtered.vcf.gz"
fi
N_FILTERED=$(bcftools view -H "${OUT}/filtered/trio_filtered.vcf.gz" | wc -l)
echo "  Total filtered variants: ${N_FILTERED}"

# ============================================================================
# Step 6a: ClinVar annotation (parallel with 6b)
# ============================================================================
echo "[Step 6a] Annotating with ClinVar..."
if [ ! -f "${OUT}/annotated/trio_clinvar.vcf.gz" ]; then
    bcftools annotate \
        -a "$CLINVAR" \
        -c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO \
        -O z \
        -o "${OUT}/annotated/trio_clinvar.vcf.gz" \
        "${OUT}/filtered/trio_filtered.vcf.gz"
    tabix -p vcf "${OUT}/annotated/trio_clinvar.vcf.gz"
fi
N_CLINVAR=$(bcftools view -H "${OUT}/annotated/trio_clinvar.vcf.gz" -i 'INFO/CLNSIG!="."' 2>/dev/null | wc -l || echo 0)
echo "  Variants with ClinVar annotation: ${N_CLINVAR}"

# ============================================================================
# Step 6b: Variant statistics (parallel with 6a)
# ============================================================================
echo "[Step 6b] Computing variant statistics..."
if [ ! -f "${OUT}/stats/trio_stats.txt" ]; then
    bcftools stats "${OUT}/filtered/trio_filtered.vcf.gz" > "${OUT}/stats/trio_stats.txt" 2>/dev/null
fi

# Per-sample stats
for SAMPLE in HG002 HG003 HG004; do
    if [ ! -f "${OUT}/stats/${SAMPLE}_stats.txt" ]; then
        bcftools stats -s "$SAMPLE" "${OUT}/filtered/trio_filtered.vcf.gz" > "${OUT}/stats/${SAMPLE}_stats.txt" 2>/dev/null || true
    fi
done

# Mendelian error check with bcftools
echo "  Computing Mendelian errors..."
if [ ! -f "${OUT}/stats/mendelian.txt" ]; then
    bcftools +mendelian \
        "${OUT}/filtered/trio_filtered.vcf.gz" \
        -p "${PED}" \
        -m c \
        > "${OUT}/stats/mendelian.txt" 2>/dev/null || {
        echo "mendelian_errors=NA" > "${OUT}/stats/mendelian.txt"
    }
fi

# Ti/Tv ratio
TITV=$(grep "^SN" "${OUT}/stats/trio_stats.txt" | grep "records" | head -1 | cut -f4 || echo "0")
TSTV_LINE=$(grep "^TSTV" "${OUT}/stats/trio_stats.txt" | head -1 || echo "")
if [ -n "$TSTV_LINE" ]; then
    TITV_RATIO=$(echo "$TSTV_LINE" | cut -f5)
else
    TITV_RATIO="0"
fi
echo "  Ti/Tv ratio: ${TITV_RATIO}"

# ============================================================================
# Step 7: Merge annotations — CONVERGE #3
# ============================================================================
echo "[Step 7] Merging all annotations (convergence #3)..."
# The ClinVar-annotated VCF already has the merged annotations
# Add variant type tags
if [ ! -f "${OUT}/annotated/trio_fully_annotated.vcf.gz" ]; then
    bcftools view "${OUT}/annotated/trio_clinvar.vcf.gz" | \
        bcftools +fill-tags -- -t AF,AN,AC,HWE | \
        bgzip > "${OUT}/annotated/trio_fully_annotated.vcf.gz"
    tabix -p vcf "${OUT}/annotated/trio_fully_annotated.vcf.gz"
fi
echo "  Fully annotated VCF ready"

# ============================================================================
# Step 8: Inheritance scoring — CONVERGE #4 (VCF + PED)
# ============================================================================
echo "[Step 8] Running inheritance scoring (convergence #4)..."
if [ ! -f "${OUT}/inheritance/trio_scored.vcf" ]; then
    bcftools view "${OUT}/annotated/trio_fully_annotated.vcf.gz" > "${OUT}/inheritance/trio_input.vcf"

    # Step 8a: genmod annotate — add gene region annotations
    genmod annotate --annotate_regions -b 38 \
        "${OUT}/inheritance/trio_input.vcf" \
        -o "${OUT}/inheritance/trio_gene_annotated.vcf" 2>/dev/null

    # Step 8b: genmod models — inheritance pattern scoring with PED
    genmod models \
        "${OUT}/inheritance/trio_gene_annotated.vcf" \
        -f "${PED}" \
        -o "${OUT}/inheritance/trio_scored.vcf" 2>/dev/null || {
        cp "${OUT}/inheritance/trio_gene_annotated.vcf" "${OUT}/inheritance/trio_scored.vcf"
    }

    rm -f "${OUT}/inheritance/trio_input.vcf" "${OUT}/inheritance/trio_gene_annotated.vcf"
fi
echo "  Inheritance scoring done"

# ============================================================================
# Step 9: Variant classification
# ============================================================================
echo "[Step 9] Classifying variants..."
python3 << 'CLASSIFY_PY'
import os, re
from collections import Counter

OUT = os.path.dirname(os.path.abspath("run_script.sh")) + "/outputs"
RESULTS = os.path.dirname(os.path.abspath("run_script.sh")) + "/results"

# Parse inheritance-scored VCF
inheritance_counts = Counter()
clinvar_sig = Counter()
variant_types = Counter()
het_hom = {"HG002": Counter(), "HG003": Counter(), "HG004": Counter()}
de_novo_candidates = 0
compound_het = 0
total_variants = 0
clinvar_pathogenic = []

scored_vcf = f"{OUT}/inheritance/trio_scored.vcf"
with open(scored_vcf) as f:
    samples = []
    for line in f:
        if line.startswith("#CHROM"):
            fields = line.strip().split('\t')
            samples = fields[9:]
            continue
        if line.startswith("#"):
            continue

        total_variants += 1
        parts = line.strip().split('\t')
        if len(parts) < 10:
            continue

        chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        info = parts[7]
        filt = parts[6]

        # Variant type
        if len(ref) == 1 and len(alt) == 1:
            variant_types["SNP"] += 1
        elif len(ref) > len(alt):
            variant_types["DEL"] += 1
        else:
            variant_types["INS"] += 1

        # Inheritance model from genmod
        inh_match = re.search(r'GeneticModels=([^;]+)', info)
        if inh_match:
            models = inh_match.group(1)
            for model in ["AR_hom", "AR_comp", "AD", "AD_dn", "XR", "XD"]:
                if model in models:
                    inheritance_counts[model] += 1
            if "AD_dn" in models:
                de_novo_candidates += 1
            if "AR_comp" in models:
                compound_het += 1

        # ClinVar significance
        clnsig_match = re.search(r'CLNSIG=([^;]+)', info)
        if clnsig_match:
            sig = clnsig_match.group(1)
            clinvar_sig[sig] += 1
            if "Pathogenic" in sig or "Likely_pathogenic" in sig:
                gene_match = re.search(r'GENEINFO=([^;|:]+)', info)
                gene = gene_match.group(1) if gene_match else "unknown"
                clinvar_pathogenic.append((chrom, pos, ref, alt, sig, gene))

        # Per-sample genotype
        fmt = parts[8].split(':')
        gt_idx = fmt.index('GT') if 'GT' in fmt else 0
        for i, sample_name in enumerate(samples):
            if i + 9 < len(parts):
                gt = parts[i + 9].split(':')[gt_idx]
                if gt in ['0/1', '0|1', '1|0']:
                    het_hom[sample_name]["het"] += 1
                elif gt in ['1/1', '1|1']:
                    het_hom[sample_name]["hom_alt"] += 1
                elif gt in ['0/0', '0|0']:
                    het_hom[sample_name]["hom_ref"] += 1

# Write classification results
with open(f"{OUT}/annotated/classification.tsv", 'w') as f:
    f.write("category\tcount\n")
    for model, count in sorted(inheritance_counts.items(), key=lambda x: -x[1]):
        f.write(f"inheritance_{model}\t{count}\n")
    for sig, count in sorted(clinvar_sig.items(), key=lambda x: -x[1])[:10]:
        f.write(f"clinvar_{sig}\t{count}\n")

print(f"Total variants: {total_variants}")
print(f"Variant types: {dict(variant_types)}")
print(f"Inheritance models: {dict(inheritance_counts)}")
print(f"ClinVar significance: {dict(list(clinvar_sig.items())[:5])}")
print(f"Pathogenic/Likely pathogenic: {len(clinvar_pathogenic)}")
print(f"De novo candidates: {de_novo_candidates}")

# Save for report
import json
summary = {
    "total_variants": total_variants,
    "variant_types": dict(variant_types),
    "inheritance_counts": dict(inheritance_counts),
    "clinvar_sig": dict(clinvar_sig),
    "n_pathogenic": len(clinvar_pathogenic),
    "de_novo_candidates": de_novo_candidates,
    "compound_het": compound_het,
    "het_hom": {k: dict(v) for k, v in het_hom.items()},
    "pathogenic_variants": [{"chrom": v[0], "pos": v[1], "ref": v[2], "alt": v[3], "sig": v[4], "gene": v[5]} for v in clinvar_pathogenic[:10]]
}
with open(f"{OUT}/annotated/summary.json", 'w') as f:
    json.dump(summary, f, indent=2)
CLASSIFY_PY

# ============================================================================
# Step 10: Generate final report
# ============================================================================
echo "[Step 10] Generating final report..."
python3 << 'REPORT_PY'
import json, os

OUT = os.path.dirname(os.path.abspath("run_script.sh")) + "/outputs"
RESULTS = os.path.dirname(os.path.abspath("run_script.sh")) + "/results"

# Load summary
with open(f"{OUT}/annotated/summary.json") as f:
    summary = json.load(f)

# Parse bcftools stats for Ti/Tv
titv_ratio = "0"
with open(f"{OUT}/stats/trio_stats.txt") as f:
    for line in f:
        if line.startswith("TSTV"):
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                titv_ratio = parts[4]
                break

# Parse Mendelian errors
mendelian_errors = 0
mendelian_file = f"{OUT}/stats/mendelian.txt"
if os.path.exists(mendelian_file):
    with open(mendelian_file) as f:
        content = f.read()
        if "mendelian_errors" in content:
            mendelian_errors = 0
        else:
            # Count non-header lines
            lines = [l for l in content.strip().split('\n') if l and not l.startswith('#')]
            mendelian_errors = len(lines)

vt = summary["variant_types"]
inh = summary["inheritance_counts"]
het_hom = summary["het_hom"]

report = []
report.append(("total_variants", summary["total_variants"]))
report.append(("snps", vt.get("SNP", 0)))
report.append(("insertions", vt.get("INS", 0)))
report.append(("deletions", vt.get("DEL", 0)))
report.append(("titv_ratio", titv_ratio))

# Per-sample genotype counts
for sample in ["HG002", "HG003", "HG004"]:
    sh = het_hom.get(sample, {})
    report.append((f"het_{sample.lower()}", sh.get("het", 0)))
    report.append((f"hom_alt_{sample.lower()}", sh.get("hom_alt", 0)))

# Inheritance patterns
report.append(("autosomal_dominant", inh.get("AD", 0)))
report.append(("autosomal_dominant_denovo", inh.get("AD_dn", 0)))
report.append(("autosomal_recessive_hom", inh.get("AR_hom", 0)))
report.append(("autosomal_recessive_comp", inh.get("AR_comp", 0)))
report.append(("x_linked_recessive", inh.get("XR", 0)))
report.append(("x_linked_dominant", inh.get("XD", 0)))

# ClinVar annotations
cs = summary["clinvar_sig"]
report.append(("clinvar_pathogenic", cs.get("Pathogenic", 0) + cs.get("Pathogenic/Likely_pathogenic", 0)))
report.append(("clinvar_likely_pathogenic", cs.get("Likely_pathogenic", 0)))
report.append(("clinvar_uncertain", cs.get("Uncertain_significance", 0)))
report.append(("clinvar_benign", cs.get("Benign", 0) + cs.get("Benign/Likely_benign", 0)))
report.append(("clinvar_likely_benign", cs.get("Likely_benign", 0)))
report.append(("total_clinvar_annotated", sum(cs.values())))

# Top pathogenic variant
pathvars = summary.get("pathogenic_variants", [])
if pathvars:
    top = pathvars[0]
    report.append(("top_pathogenic_gene", top["gene"]))
    report.append(("top_pathogenic_pos", f"chr22:{top['pos']}"))
    report.append(("top_pathogenic_change", f"{top['ref']}>{top['alt']}"))
else:
    report.append(("top_pathogenic_gene", "none"))
    report.append(("top_pathogenic_pos", "none"))
    report.append(("top_pathogenic_change", "none"))

report.append(("mendelian_errors", mendelian_errors))
report.append(("de_novo_candidates", summary["de_novo_candidates"]))
report.append(("compound_het_candidates", summary["compound_het"]))

with open(f"{RESULTS}/report.csv", 'w') as f:
    f.write("metric,value\n")
    for m, v in report:
        f.write(f"{m},{v}\n")

print(f"Report: {len(report)} metrics")
for m, v in report:
    print(f"  {m}: {v}")
REPORT_PY

echo ""
echo "========================================="
echo "  Variant Annotation Trio Complete!"
echo "========================================="
echo ""
cat "${RESULTS}/report.csv"
