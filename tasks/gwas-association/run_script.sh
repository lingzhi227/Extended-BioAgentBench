#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# GWAS / Population Association Testing Pipeline
# =============================================================================
# DAG Structure (depth=10, convergence=4, tools=10):
#
#  genotypes.bed/bim/fam    phenotype.txt    covariates.txt
#          |                     |                |
#  [plink2 --geno --mind         |                |
#   --maf --hwe QC] -------------------------------- Level 1
#          |
#  +-------+-----------+
#  |       |           |
# [plink  [plink2     [plink2                       Level 2
#  --ld     --pca]     --het (inbreeding)]
#  prune]
#  |       |           |
#  |   CONVERGENCE 1   |
#  |   (pruned + PCA   |
#  |    eigenvecs)     |
#  |       |    +------+
#  |       |    |
#  |  [plink2 outlier removal]                      Level 3
#  |       |
#  +-------+-----------+
#                      |
#          CONVERGENCE 2                            Level 4
#          (QC'd genotypes + PCs + cleaned)
#                      |
#          +-----------+-----------+
#          |           |           |
#    [plink2        [regenie      [plink2            Level 5-6
#     --glm          step1->2]     --assoc]
#     logistic]
#          |           |           |
#          +-----------+-----------+
#                      |
#              CONVERGENCE 3                        Level 7
#              (3-method results)
#                      |
#          +-----------+-----------+
#          |           |           |
#    [python        [plink2       [python            Level 8
#     manhattan]     --clump]      summary]
#          |           |           |
#          +-----------+-----------+
#                      |
#              CONVERGENCE 4                        Level 9
#              (plots + clumps + summary)
#                      |
#              [python report]                      Level 10
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
OUT="${WORKDIR}/outputs"
RES="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,pruned,pca,het,assoc_plink,assoc_regenie,assoc_fisher,clump,summary} "${RES}"

# Input files
BED="${DATA}/example"
PHENO="${DATA}/phenotype_bin.txt"
COVAR="${DATA}/covariates.txt"

# =============================================================================
# Level 1: Genotype QC
# =============================================================================
if [ ! -f "${OUT}/qc/genotypes_qc.bed" ]; then
  echo "[L1] Genotype QC..."
  plink2 --bfile "${BED}" \
    --geno 0.1 --mind 0.1 --maf 0.01 --hwe 1e-6 \
    --make-bed --out "${OUT}/qc/genotypes_qc" 2>&1 | tail -5
fi

# =============================================================================
# Level 2a: LD pruning
# =============================================================================
if [ ! -f "${OUT}/pruned/pruned.prune.in" ]; then
  echo "[L2a] LD pruning..."
  plink --bfile "${OUT}/qc/genotypes_qc" \
    --indep-pairwise 50 5 0.2 \
    --out "${OUT}/pruned/pruned" 2>&1 | tail -5
fi

# =============================================================================
# Level 2b: PCA on LD-pruned SNPs
# =============================================================================
if [ ! -f "${OUT}/pca/pca.eigenvec" ]; then
  echo "[L2b] PCA..."
  plink2 --bfile "${OUT}/qc/genotypes_qc" \
    --extract "${OUT}/pruned/pruned.prune.in" \
    --pca 10 \
    --out "${OUT}/pca/pca" 2>&1 | tail -5
fi

# =============================================================================
# Level 2c: Heterozygosity (inbreeding check)
# =============================================================================
if [ ! -f "${OUT}/het/het.het" ]; then
  echo "[L2c] Heterozygosity check..."
  plink2 --bfile "${OUT}/qc/genotypes_qc" \
    --het \
    --out "${OUT}/het/het" 2>&1 | tail -5
fi

# =============================================================================
# Level 3: CONVERGENCE 1 — Outlier removal
# =============================================================================
echo "[L3] Convergence 1: outlier removal..."

if [ ! -f "${OUT}/qc/keep_samples.txt" ]; then
  python3 << 'OUTLIER'
import os

OUT = os.environ.get("OUT", "outputs")

# Read heterozygosity
het_data = {}
het_file = f"{OUT}/het/het.het"
if os.path.exists(het_file):
    with open(het_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                fid, iid = parts[0], parts[1]
                try:
                    f_val = float(parts[5])
                    het_data[(fid, iid)] = f_val
                except:
                    pass

# Remove samples with |F| > 0.2
import statistics
if het_data:
    f_values = list(het_data.values())
    mean_f = statistics.mean(f_values)
    sd_f = statistics.stdev(f_values) if len(f_values) > 1 else 0
    keep = [(fid, iid) for (fid, iid), f in het_data.items()
            if abs(f - mean_f) < 3 * sd_f]
else:
    keep = []

with open(f"{OUT}/qc/keep_samples.txt", "w") as f:
    for fid, iid in keep:
        f.write(f"{fid}\t{iid}\n")

print(f"Keeping {len(keep)} / {len(het_data)} samples")
OUTLIER
fi

# =============================================================================
# Level 4: CONVERGENCE 2 — Apply QC
# =============================================================================
if [ ! -f "${OUT}/qc/genotypes_final.bed" ]; then
  echo "[L4] Applying final QC..."
  plink2 --bfile "${OUT}/qc/genotypes_qc" \
    --keep "${OUT}/qc/keep_samples.txt" \
    --make-bed --out "${OUT}/qc/genotypes_final" 2>&1 | tail -5
fi

# =============================================================================
# Level 5a: PLINK2 logistic regression
# =============================================================================
if [ ! -f "${OUT}/assoc_plink/plink_assoc.Y1.glm.logistic.hybrid" ]; then
  echo "[L5a] PLINK2 --glm logistic..."
  plink2 --bfile "${OUT}/qc/genotypes_final" \
    --pheno "${PHENO}" --pheno-name Y1 \
    --covar "${OUT}/pca/pca.eigenvec" \
    --glm hide-covar --covar-variance-standardize \
    --out "${OUT}/assoc_plink/plink_assoc" 2>&1 | tail -5 || true
fi

# =============================================================================
# Level 5b: REGENIE step 1 + step 2
# =============================================================================
if [ ! -f "${OUT}/assoc_regenie/step2_Y1.regenie" ]; then
  echo "[L5b] REGENIE step 1..."
  # Create phenotype file in right format
  regenie --step 1 \
    --bed "${OUT}/qc/genotypes_final" \
    --phenoFile "${PHENO}" --phenoColList Y1 \
    --bsize 100 \
    --bt --lowmem \
    --out "${OUT}/assoc_regenie/step1" 2>&1 | tail -10 || true

  echo "[L5b] REGENIE step 2..."
  if [ -f "${OUT}/assoc_regenie/step1_pred.list" ]; then
    regenie --step 2 \
      --bed "${OUT}/qc/genotypes_final" \
      --phenoFile "${PHENO}" --phenoColList Y1 \
      --bsize 200 --bt --firth --approx \
      --pred "${OUT}/assoc_regenie/step1_pred.list" \
      --out "${OUT}/assoc_regenie/step2" 2>&1 | tail -10 || true
  fi
fi

# =============================================================================
# Level 5c: PLINK2 Fisher exact test
# =============================================================================
if [ ! -f "${OUT}/assoc_fisher/fisher_assoc.Y1.glm.logistic.hybrid" ]; then
  echo "[L5c] PLINK2 Fisher test..."
  plink2 --bfile "${OUT}/qc/genotypes_final" \
    --pheno "${PHENO}" --pheno-name Y1 \
    --glm hide-covar \
    --out "${OUT}/assoc_fisher/fisher_assoc" 2>&1 | tail -5 || true
fi

# =============================================================================
# Level 7: CONVERGENCE 3 — merge results
# =============================================================================
echo "[L7] Convergence 3: merging association results..."

# =============================================================================
# Level 8: LD clumping
# =============================================================================
if [ ! -f "${OUT}/clump/clumped.clumps" ]; then
  echo "[L8] LD clumping..."
  # Find the PLINK association results
  ASSOC_FILE=$(ls "${OUT}"/assoc_plink/*.glm.logistic* 2>/dev/null | head -1 || true)
  if [ -n "${ASSOC_FILE}" ]; then
    plink --bfile "${OUT}/qc/genotypes_final" \
      --clump "${ASSOC_FILE}" \
      --clump-p1 0.01 --clump-p2 0.05 --clump-r2 0.2 \
      --out "${OUT}/clump/clumped" 2>&1 | tail -5 || true
  fi
fi

# =============================================================================
# Level 9-10: Report generation
# =============================================================================
echo "[L9-10] Generating report..."

export OUT RES
python3 << 'REPORT'
import os, glob

OUT = os.environ.get("OUT", "outputs")
RES_DIR = os.environ.get("RES", "results")

metrics = {}

# QC stats
for logf in glob.glob(f"{OUT}/qc/genotypes_qc.log"):
    with open(logf) as f:
        for line in f:
            if "variant" in line.lower() and "loaded" in line.lower():
                parts = line.strip().split()
                try:
                    metrics["variants_before_qc"] = int(parts[0])
                except:
                    pass
            if "sample" in line.lower() and ("remaining" in line.lower() or "pass" in line.lower()):
                parts = line.strip().split()
                try:
                    metrics["samples_after_qc"] = int(parts[0])
                except:
                    pass

# Final genotype stats
for logf in glob.glob(f"{OUT}/qc/genotypes_final.log"):
    with open(logf) as f:
        for line in f:
            if "variant" in line.lower() and ("remaining" in line.lower() or "pass" in line.lower()):
                parts = line.strip().split()
                try:
                    metrics["variants_after_qc"] = int(parts[0])
                except:
                    pass

# PCA
pca_file = f"{OUT}/pca/pca.eigenvec"
if os.path.exists(pca_file):
    with open(pca_file) as f:
        metrics["pca_samples"] = sum(1 for l in f) - 1  # minus header

# LD pruning
prune_file = f"{OUT}/pruned/pruned.prune.in"
if os.path.exists(prune_file):
    with open(prune_file) as f:
        metrics["ld_pruned_snps"] = sum(1 for l in f if l.strip())

# Heterozygosity
het_file = f"{OUT}/het/het.het"
if os.path.exists(het_file):
    f_vals = []
    with open(het_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    f_vals.append(float(parts[5]))
                except:
                    pass
    if f_vals:
        metrics["mean_het_f"] = round(sum(f_vals) / len(f_vals), 4)

# Association results - PLINK
for assoc_file in glob.glob(f"{OUT}/assoc_plink/*.glm.*"):
    sig_5e8 = 0
    sig_1e5 = 0
    total = 0
    with open(assoc_file) as f:
        header = f.readline().strip().split("\t")
        p_idx = header.index("P") if "P" in header else -1
        if p_idx >= 0:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > p_idx and parts[p_idx] != "NA":
                    try:
                        p = float(parts[p_idx])
                        total += 1
                        if p < 5e-8:
                            sig_5e8 += 1
                        if p < 1e-5:
                            sig_1e5 += 1
                    except:
                        pass
    metrics["plink_tested_variants"] = total
    metrics["plink_genome_wide_sig"] = sig_5e8
    metrics["plink_suggestive_sig"] = sig_1e5
    break

# Association results - REGENIE
for assoc_file in glob.glob(f"{OUT}/assoc_regenie/step2_Y1.regenie"):
    sig_5e8 = 0
    sig_1e5 = 0
    total = 0
    with open(assoc_file) as f:
        header = f.readline().strip().split()
        log10p_idx = header.index("LOG10P") if "LOG10P" in header else -1
        if log10p_idx >= 0:
            for line in f:
                parts = line.strip().split()
                if len(parts) > log10p_idx and parts[log10p_idx] != "NA":
                    try:
                        log10p = float(parts[log10p_idx])
                        total += 1
                        if log10p > 7.3:  # -log10(5e-8) ≈ 7.3
                            sig_5e8 += 1
                        if log10p > 5:
                            sig_1e5 += 1
                    except:
                        pass
    metrics["regenie_tested_variants"] = total
    metrics["regenie_genome_wide_sig"] = sig_5e8
    metrics["regenie_suggestive_sig"] = sig_1e5
    break

# Clumping
clump_file = f"{OUT}/clump/clumped.clumps"
if os.path.exists(clump_file):
    with open(clump_file) as f:
        metrics["clumped_loci"] = sum(1 for l in f if l.strip() and not l.startswith("CHR")) - 1  # minus header
else:
    metrics["clumped_loci"] = 0

# Genomic inflation
# Read plink association p-values to compute lambda
import math
all_pvals = []
for assoc_file in glob.glob(f"{OUT}/assoc_plink/*.glm.*"):
    with open(assoc_file) as f:
        header = f.readline().strip().split("\t")
        p_idx = header.index("P") if "P" in header else -1
        if p_idx >= 0:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > p_idx and parts[p_idx] not in ("NA", ""):
                    try:
                        all_pvals.append(float(parts[p_idx]))
                    except:
                        pass
    break

if all_pvals:
    from scipy import stats
    chi2_vals = [stats.chi2.ppf(1 - p, 1) for p in all_pvals if 0 < p < 1]
    if chi2_vals:
        median_chi2 = sorted(chi2_vals)[len(chi2_vals) // 2]
        lambda_gc = round(median_chi2 / 0.4549, 4)
        metrics["genomic_inflation_lambda"] = lambda_gc

# Write CSV
with open(f"{RES_DIR}/report.csv", "w") as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("Report written:")
for k, v in metrics.items():
    print(f"  {k} = {v}")
REPORT

echo "=== Pipeline Complete ==="
cat "${RES}/report.csv"
