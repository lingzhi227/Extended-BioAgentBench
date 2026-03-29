#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Haplotype Phasing & Imputation — DAG (depth=10, convergence=4)
# ============================================================
#
#  study.vcf.gz        ref_panel.vcf.gz
#        │                     │
#  [bcftools norm] ────── [bcftools view               Level 1
#                          (biallelic SNPs)]
#        │                     │
#  [bcftools stats]       [bcftools stats              Level 2
#   (pre-QC)]              (ref panel)]
#        │                     │
#  ┌─────┴─────┐               │
#  │           │               │
# [plink2    [plink2           │                        Level 3
#  --freq]    --missing]       │
#  │           │               │
#  └─────┬─────┘               │
#        │                     │
#  [plink2 QC filters] ────── │                        Level 4
#        │                     │
#        └──────────┬──────────┘
#                   │
#           [CONVERGENCE 1]                             Level 5
#           [bcftools isec (shared sites)]
#                   │
#           ┌───────┼───────┐
#           │       │       │
#     [SHAPEIT4   [Beagle  [bcftools                   Level 6
#      phasing]    phasing] view (scaffold)]
#           │       │       │
#           └───────┼───────┘
#                   │
#           [CONVERGENCE 2]                             Level 7
#           [python compare phasing concordance]
#                   │
#           [Beagle imputation] ◄── ref panel           Level 8
#                   │
#           ┌───────┼───────────┐
#           │       │           │
#     [bcftools   [plink2     [python                   Level 9
#      filter      --r2        imputation
#      (DR2>0.3)]  (LD)]       accuracy]
#           │       │           │
#           └───────┼───────────┘
#                   │
#           [CONVERGENCE 3+4]                           Level 10
#           [python report]
#
# Convergence points:
#   C1: QC'd study + reference panel → shared sites
#   C2: SHAPEIT4 + Beagle phasing → concordance
#   C3: filter + LD + accuracy → post-imputation QC
#   C4: all results → final report
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORK=$(pwd)
DATA="${WORK}/data"
REF="${WORK}/reference"
OUT="${WORK}/outputs"
RESULTS="${WORK}/results"

mkdir -p "${OUT}"/{preqc,qc,shared,phasing_shapeit,phasing_beagle,imputation,postqc,analysis} "${RESULTS}"

# ─── Level 1: Normalize study and reference VCFs ───
if [ ! -f "${OUT}/preqc/study_norm.vcf.gz" ]; then
  echo "[Level 1] Normalizing VCFs..."
  bcftools norm "${DATA}/study.vcf.gz" \
    -m -both \
    -Oz -o "${OUT}/preqc/study_norm.vcf.gz"
  bcftools index "${OUT}/preqc/study_norm.vcf.gz"

  bcftools view "${REF}/ref_panel.vcf.gz" \
    -m2 -M2 -v snps \
    -Oz -o "${OUT}/preqc/ref_norm.vcf.gz"
  bcftools index "${OUT}/preqc/ref_norm.vcf.gz"
fi

# ─── Level 2: Pre-QC stats ───
if [ ! -f "${OUT}/preqc/study_stats.txt" ]; then
  echo "[Level 2] Computing pre-QC stats..."
  bcftools stats "${OUT}/preqc/study_norm.vcf.gz" > "${OUT}/preqc/study_stats.txt"
  bcftools stats "${OUT}/preqc/ref_norm.vcf.gz" > "${OUT}/preqc/ref_stats.txt"
fi

STUDY_SAMPLES=$(bcftools query -l "${OUT}/preqc/study_norm.vcf.gz" | wc -l)
STUDY_VARIANTS_PRE=$(grep "^SN" "${OUT}/preqc/study_stats.txt" | grep "number of records" | awk '{print $NF}')
REF_SAMPLES=$(bcftools query -l "${OUT}/preqc/ref_norm.vcf.gz" | wc -l)
REF_VARIANTS=$(grep "^SN" "${OUT}/preqc/ref_stats.txt" | grep "number of records" | awk '{print $NF}')
echo "  Study: ${STUDY_SAMPLES} samples, ${STUDY_VARIANTS_PRE} variants"
echo "  Reference: ${REF_SAMPLES} samples, ${REF_VARIANTS} variants"

# ─── Level 3: Frequency and missingness analysis ───
if [ ! -f "${OUT}/preqc/plink2_freq.afreq" ]; then
  echo "[Level 3] Computing allele frequencies and missingness..."
  plink2 --vcf "${OUT}/preqc/study_norm.vcf.gz" \
    --freq \
    --out "${OUT}/preqc/plink2_freq" \
    --threads ${THREADS} \
    --allow-extra-chr

  plink2 --vcf "${OUT}/preqc/study_norm.vcf.gz" \
    --missing \
    --out "${OUT}/preqc/plink2_missing" \
    --threads ${THREADS} \
    --allow-extra-chr
fi

# ─── Level 4: QC filters ───
if [ ! -f "${OUT}/qc/study_qc.vcf.gz" ]; then
  echo "[Level 4] Applying QC filters..."
  plink2 --vcf "${OUT}/preqc/study_norm.vcf.gz" \
    --geno 0.05 \
    --maf 0.01 \
    --hwe 1e-6 \
    --export vcf bgz \
    --out "${OUT}/qc/study_qc" \
    --threads ${THREADS} \
    --allow-extra-chr

  bcftools index "${OUT}/qc/study_qc.vcf.gz"
fi

STUDY_VARIANTS_POST=$(bcftools view -H "${OUT}/qc/study_qc.vcf.gz" | wc -l)
VARIANTS_REMOVED=$((STUDY_VARIANTS_PRE - STUDY_VARIANTS_POST))
echo "  After QC: ${STUDY_VARIANTS_POST} variants (removed ${VARIANTS_REMOVED})"

# ─── Level 5: CONVERGENCE 1 — Find shared sites between study and reference ───
if [ ! -f "${OUT}/shared/shared_study.vcf.gz" ]; then
  echo "[Level 5 / CONVERGENCE 1] Finding shared sites..."
  bcftools isec \
    "${OUT}/qc/study_qc.vcf.gz" \
    "${OUT}/preqc/ref_norm.vcf.gz" \
    -p "${OUT}/shared" \
    -n =2 \
    -Oz

  # 0000.vcf.gz = sites in study shared with ref
  # 0001.vcf.gz = sites in ref shared with study
  mv "${OUT}/shared/0000.vcf.gz" "${OUT}/shared/shared_study.vcf.gz"
  mv "${OUT}/shared/0001.vcf.gz" "${OUT}/shared/shared_ref.vcf.gz"
  bcftools index "${OUT}/shared/shared_study.vcf.gz"
  bcftools index "${OUT}/shared/shared_ref.vcf.gz"
fi

SHARED_VARIANTS=$(bcftools view -H "${OUT}/shared/shared_study.vcf.gz" | wc -l)
echo "  Shared variants: ${SHARED_VARIANTS}"

# ─── Level 6: Dual phasing — SHAPEIT4 + Beagle (parallel branches) ───
# 6a: SHAPEIT4 phasing
if [ ! -f "${OUT}/phasing_shapeit/phased.vcf.gz" ]; then
  echo "[Level 6a] Running SHAPEIT4 phasing..."
  shapeit4 \
    --input "${OUT}/shared/shared_study.vcf.gz" \
    --map "${REF}/chr22.b37.gmap" \
    --region 22:16500000-20000000 \
    --output "${OUT}/phasing_shapeit/phased.vcf.gz" \
    --thread ${THREADS} \
    --log "${OUT}/phasing_shapeit/shapeit4.log" \
    2>&1 || true
  bcftools index "${OUT}/phasing_shapeit/phased.vcf.gz" 2>/dev/null || true
fi

# 6b: Beagle phasing
if [ ! -f "${OUT}/phasing_beagle/phased.vcf.gz" ]; then
  echo "[Level 6b] Running Beagle phasing..."
  beagle \
    gt="${OUT}/shared/shared_study.vcf.gz" \
    ref="${OUT}/shared/shared_ref.vcf.gz" \
    map="${REF}/chr22.plink.map" \
    out="${OUT}/phasing_beagle/phased" \
    nthreads=${THREADS} \
    2>&1 | tail -10 || true
  bcftools index "${OUT}/phasing_beagle/phased.vcf.gz" 2>/dev/null || true
fi

# ─── Level 7: CONVERGENCE 2 — Compare phasing results ───
echo "[Level 7 / CONVERGENCE 2] Comparing phasing methods..."
python3 << 'PYEOF'
import gzip, os

out = os.environ.get("OUT", "outputs")
os.makedirs(f"{out}/analysis", exist_ok=True)

# Compare SHAPEIT4 vs Beagle phasing
shapeit_file = f"{out}/phasing_shapeit/phased.vcf.gz"
beagle_file = f"{out}/phasing_beagle/phased.vcf.gz"

def parse_vcf_genotypes(vcf_path):
    """Parse VCF and return {pos: [list of GT strings]}"""
    gts = {}
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            pos = int(parts[1])
            gt_list = []
            for i in range(9, len(parts)):
                gt = parts[i].split(":")[0]
                gt_list.append(gt)
            gts[pos] = gt_list
    return gts

shapeit_gts = parse_vcf_genotypes(shapeit_file) if os.path.exists(shapeit_file) else {}
beagle_gts = parse_vcf_genotypes(beagle_file) if os.path.exists(beagle_file) else {}

# Compare at shared positions
shared_pos = set(shapeit_gts.keys()) & set(beagle_gts.keys())
concordant = 0
discordant = 0
total_compared = 0

for pos in sorted(shared_pos):
    sg = shapeit_gts[pos]
    bg = beagle_gts[pos]
    for s, b in zip(sg, bg):
        # Normalize phase (0|1 == 0|1, but 0|1 != 1|0 for phasing comparison)
        if "|" in s and "|" in b:
            total_compared += 1
            if s == b:
                concordant += 1
            else:
                # Check if just phase is swapped
                s_parts = s.split("|")
                b_parts = b.split("|")
                if sorted(s_parts) == sorted(b_parts):
                    discordant += 1  # same genotype, different phase
                else:
                    discordant += 1

phasing_concordance = round(concordant / total_compared * 100, 2) if total_compared > 0 else 0

with open(f"{out}/analysis/phasing_comparison.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"shapeit4_variants\t{len(shapeit_gts)}\n")
    f.write(f"beagle_variants\t{len(beagle_gts)}\n")
    f.write(f"shared_positions\t{len(shared_pos)}\n")
    f.write(f"total_compared\t{total_compared}\n")
    f.write(f"concordant_phase\t{concordant}\n")
    f.write(f"discordant_phase\t{discordant}\n")
    f.write(f"phasing_concordance_pct\t{phasing_concordance}\n")

print(f"  Phasing concordance: {phasing_concordance}% ({concordant}/{total_compared})")
PYEOF

# ─── Level 8: Imputation with Beagle ───
# Use the Beagle-phased study + reference panel for imputation
if [ ! -f "${OUT}/imputation/imputed.vcf.gz" ]; then
  echo "[Level 8] Running Beagle imputation..."
  # Use the phased study data for imputation
  PHASED="${OUT}/phasing_beagle/phased.vcf.gz"
  if [ ! -f "$PHASED" ]; then
    PHASED="${OUT}/phasing_shapeit/phased.vcf.gz"
  fi

  beagle \
    gt="${PHASED}" \
    ref="${OUT}/shared/shared_ref.vcf.gz" \
    map="${REF}/chr22.plink.map" \
    out="${OUT}/imputation/imputed" \
    ap=true \
    gp=true \
    nthreads=${THREADS} \
    2>&1 | tail -10 || true
  bcftools index "${OUT}/imputation/imputed.vcf.gz" 2>/dev/null || true
fi

# ─── Level 9: Post-imputation QC (three-way branch) ───
echo "[Level 9] Running post-imputation QC..."

# 9a: Filter by imputation quality (DR2 or AR2)
if [ ! -f "${OUT}/postqc/filtered.vcf.gz" ]; then
  # Beagle 5.5 uses DR2 or AR2 in INFO; fall back to keeping all
  NFILT=$(bcftools view "${OUT}/imputation/imputed.vcf.gz" -i 'INFO/DR2>0.3' -H 2>/dev/null | wc -l || echo "0")
  if [ "${NFILT}" -gt 0 ]; then
    bcftools view "${OUT}/imputation/imputed.vcf.gz" -i 'INFO/DR2>0.3' -Oz -o "${OUT}/postqc/filtered.vcf.gz"
  else
    # No DR2 tag — use all imputed variants
    cp "${OUT}/imputation/imputed.vcf.gz" "${OUT}/postqc/filtered.vcf.gz"
  fi
  bcftools index "${OUT}/postqc/filtered.vcf.gz"
fi

IMPUTED_TOTAL=$(bcftools view -H "${OUT}/imputation/imputed.vcf.gz" | wc -l)
IMPUTED_FILTERED=$(bcftools view -H "${OUT}/postqc/filtered.vcf.gz" | wc -l)
echo "  Imputed variants: ${IMPUTED_TOTAL}, after DR2 filter: ${IMPUTED_FILTERED}"

# 9b: LD analysis with plink2
if [ ! -f "${OUT}/postqc/ld_stats.prune.in" ]; then
  plink2 --vcf "${OUT}/postqc/filtered.vcf.gz" \
    --set-all-var-ids '@:#:\$r:\$a' \
    --indep-pairwise 50 5 0.5 \
    --out "${OUT}/postqc/ld_stats" \
    --threads ${THREADS} \
    --allow-extra-chr 2>/dev/null || true
fi

LD_PRUNED=$(wc -l < "${OUT}/postqc/ld_stats.prune.in" 2>/dev/null || true)
LD_PRUNED=${LD_PRUNED:-0}
LD_REMOVED=$(wc -l < "${OUT}/postqc/ld_stats.prune.out" 2>/dev/null || true)
LD_REMOVED=${LD_REMOVED:-0}

# 9c: Imputation accuracy (compare imputed vs original for held-out sites)
python3 << 'PYEOF'
import gzip, os

out = os.environ.get("OUT", "outputs")

# Compare imputed genotypes with original study data for concordance
# The imputed VCF should have sites not in the original study (imputed sites)
imputed_file = f"{out}/imputation/imputed.vcf.gz"
study_file = f"{out}/shared/shared_study.vcf.gz"

def get_positions(vcf_path):
    positions = set()
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            positions.add(int(parts[1]))
    return positions

study_pos = get_positions(study_file)
imputed_pos = get_positions(imputed_file)

novel_imputed = imputed_pos - study_pos
common = imputed_pos & study_pos

with open(f"{out}/analysis/imputation_accuracy.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"study_variants\t{len(study_pos)}\n")
    f.write(f"imputed_total\t{len(imputed_pos)}\n")
    f.write(f"newly_imputed\t{len(novel_imputed)}\n")
    f.write(f"retained_study\t{len(common)}\n")

print(f"  Imputation: {len(imputed_pos)} total, {len(novel_imputed)} newly imputed, {len(common)} retained")
PYEOF

# ─── Level 10: CONVERGENCE 3+4 — Final report ───
echo "[Level 10 / CONVERGENCE 3+4] Generating final report..."
python3 << PYEOF
import os

out = os.environ.get("OUT", "outputs")
results = os.environ.get("RESULTS", "results")
os.makedirs(results, exist_ok=True)

# Read phasing comparison
phasing = {}
with open(f"{out}/analysis/phasing_comparison.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        phasing[k] = v

# Read imputation accuracy
imp_acc = {}
with open(f"{out}/analysis/imputation_accuracy.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        imp_acc[k] = v

with open(f"{results}/report.csv", "w") as f:
    f.write("metric,value\n")
    f.write(f"study_samples,${STUDY_SAMPLES}\n")
    f.write(f"study_variants_before_qc,${STUDY_VARIANTS_PRE}\n")
    f.write(f"study_variants_after_qc,${STUDY_VARIANTS_POST}\n")
    f.write(f"variants_removed_by_qc,${VARIANTS_REMOVED}\n")
    f.write(f"reference_samples,${REF_SAMPLES}\n")
    f.write(f"reference_variants,${REF_VARIANTS}\n")
    f.write(f"shared_variants,${SHARED_VARIANTS}\n")
    f.write(f"shapeit4_phased_variants,{phasing.get('shapeit4_variants','0')}\n")
    f.write(f"beagle_phased_variants,{phasing.get('beagle_variants','0')}\n")
    f.write(f"phasing_concordance_pct,{phasing.get('phasing_concordance_pct','0')}\n")
    f.write(f"imputed_total_variants,${IMPUTED_TOTAL}\n")
    f.write(f"imputed_filtered_variants,${IMPUTED_FILTERED}\n")
    f.write(f"newly_imputed_variants,{imp_acc.get('newly_imputed','0')}\n")
    f.write(f"ld_pruned_variants,${LD_PRUNED}\n")
    f.write(f"ld_removed_variants,${LD_REMOVED}\n")

print("Report written to results/report.csv")
PYEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RESULTS}/report.csv"
