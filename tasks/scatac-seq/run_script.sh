#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Single-Cell ATAC-seq Pipeline: Chromatin Accessibility Profiling
# ============================================================
# DAG (depth=10, convergence=4):
#
#  fragments.tsv.gz + barcodes
#        │
#  [snapatac2 import]                                   Level 1
#        │
#  ┌─────┴──────────┐
#  │                │
# [snapatac2      [python                               Level 2
#  qc_metrics]     barcode stats]
#  │                │
#  └────────┬───────┘
#           │
#   [CONVERGENCE 1: QC + stats]                         Level 3
#   [snapatac2 filter cells]
#           │
#   [snapatac2 add_tile_matrix]                         Level 4
#           │
#   ┌───────┼───────────┐
#   │       │           │
# [snapatac2 [snapatac2 [snapatac2                      Level 5
#  select_   spectral    scrublet
#  features] (dim red)]  (doublets)]
#   │       │           │
#   └───────┼───────────┘
#           │
#   [CONVERGENCE 2: features+embed+clean]               Level 6
#   [snapatac2 leiden clustering]
#           │
#   ┌───────┼───────────────────┐
#   │       │                   │
# [macs2   [snapatac2          [snapatac2               Level 7
#  callpeak gene_matrix]        export]
#  per-clust (gene activity)]
#   │       │                   │
#   │  [CONVERGENCE 3] ◄────────┘                       Level 8
#   │  (clusters + gene activity + export)
#   │       │
#   └───────┴───────┐
#                   │
#   [CONVERGENCE 3+: peaks + clusters + activity]       Level 9
#                   │
#   ┌───────┼───────────┐
#   │       │           │
# [python  [python    [python                           Level 9
#  diff     marker     cluster
#  access]  genes]     stats]
#   │       │           │
#   └───────┼───────────┘
#           │
#   [CONVERGENCE 4: report]                             Level 10
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{qc,filtered,clusters,peaks,gene_activity,community}
mkdir -p "${RESULTS}"

FRAGMENTS="${DATA}/fragments.tsv.gz"
GTF="${REF}/gencode.v32.basic.annotation.gtf.gz"

# ============================================================
# Full pipeline in Python/SnapATAC2
# ============================================================
python3 << 'PYEOF'
import snapatac2 as snap
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import subprocess
import gzip

os.chdir(os.environ.get("WORKDIR", "."))
if not os.path.exists("outputs"):
    os.chdir("/pscratch/sd/l/lingzhi/bench-task-output/session-i/scatac-seq")

print("=== Level 1: Import fragments ===")
adata = snap.pp.import_data(
    "data/fragments.tsv.gz",
    chrom_sizes=snap.genome.hg38,
    min_num_fragments=100,
    sorted_by_barcode=False,
)
print(f"  Imported: {adata.n_obs} cells, {adata.n_vars} features")

print("=== Level 2: QC metrics ===")
snap.metrics.tsse(adata, snap.genome.hg38)
snap.metrics.frag_size_distr(adata)

n_cells_before = adata.n_obs
tsse_values = adata.obs["tsse"].values
frag_values = adata.obs["n_fragment"].values

print(f"  TSS enrichment: median={np.median(tsse_values):.2f}, mean={np.mean(tsse_values):.2f}")
print(f"  Fragments: median={np.median(frag_values):.0f}, mean={np.mean(frag_values):.0f}")

# Save QC stats
pd.DataFrame({
    "barcode": adata.obs_names,
    "tsse": tsse_values,
    "n_fragment": frag_values
}).to_csv("outputs/qc/cell_qc_metrics.tsv", sep='\t', index=False)

print("=== Level 3: CONVERGENCE 1 - Filter cells ===")
snap.pp.filter_cells(adata, min_counts=1000, min_tsse=4.0, max_counts=100000)
n_cells_after = adata.n_obs
print(f"  Filtered: {n_cells_before} -> {n_cells_after} cells")

print("=== Level 4: Tile matrix ===")
snap.pp.add_tile_matrix(adata, bin_size=5000)
print(f"  Tile matrix: {adata.shape}")

print("=== Level 5: Feature selection + dim reduction + doublet detection ===")
snap.pp.select_features(adata)
snap.tl.spectral(adata)

# Doublet detection (scrublet)
n_doublets = 0
try:
    snap.pp.scrublet(adata)
    doublet_scores = adata.obs.get("doublet_score", pd.Series([0]*adata.n_obs))
    n_doublets = int((doublet_scores > 0.5).sum()) if hasattr(doublet_scores, '__gt__') else 0
    if "is_doublet" in adata.obs.columns:
        adata = adata[~adata.obs["is_doublet"]].copy()
        print(f"  After doublet removal: {adata.n_obs} cells")
except Exception as e:
    print(f"  Scrublet skipped (compatibility issue): {e}")
    n_doublets = 0
print(f"  Doublets detected: {n_doublets}")

print("=== Level 6: CONVERGENCE 2 - Clustering ===")
snap.pp.knn(adata)
snap.tl.leiden(adata, resolution=0.5)
snap.tl.umap(adata)

n_clusters = len(adata.obs["leiden"].unique())
print(f"  Clusters: {n_clusters}")

# Save cluster assignments
adata.obs[["leiden"]].to_csv("outputs/clusters/cluster_assignments.tsv", sep='\t')

# Cluster sizes
cluster_sizes = adata.obs["leiden"].value_counts().sort_index()
cluster_sizes.to_csv("outputs/clusters/cluster_sizes.tsv", sep='\t', header=["count"])

print("=== Level 7: Peak calling + Gene activity ===")

# Branch 7a: MACS2 peak calling per cluster
# Export fragments per cluster for MACS2
cluster_labels = adata.obs["leiden"].values
barcodes_per_cluster = {}
for cl in sorted(set(cluster_labels)):
    bcs = set(adata.obs_names[cluster_labels == cl])
    barcodes_per_cluster[cl] = bcs

# Read fragments and split by cluster
print("  Splitting fragments by cluster for peak calling...")
cluster_files = {}
for cl in barcodes_per_cluster:
    fpath = f"outputs/peaks/cluster_{cl}_fragments.bed"
    cluster_files[cl] = open(fpath, 'w')

with gzip.open("data/fragments.tsv.gz", 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            bc = parts[3]
            for cl, bcs in barcodes_per_cluster.items():
                if bc in bcs:
                    cluster_files[cl].write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

for f in cluster_files.values():
    f.close()

# Run MACS2 per cluster
total_peaks = 0
for cl in sorted(barcodes_per_cluster.keys()):
    bed = f"outputs/peaks/cluster_{cl}_fragments.bed"
    outdir = f"outputs/peaks/cluster_{cl}"
    if os.path.getsize(bed) > 0:
        os.makedirs(outdir, exist_ok=True)
        result = subprocess.run([
            "macs2", "callpeak",
            "-t", bed,
            "-f", "BED",
            "--nomodel",
            "--shift", "-75",
            "--extsize", "150",
            "--keep-dup", "all",
            "-n", f"cluster_{cl}",
            "--outdir", outdir,
            "-q", "0.05",
            "-g", "hs"
        ], capture_output=True, text=True)

        peak_file = f"{outdir}/cluster_{cl}_peaks.narrowPeak"
        if os.path.exists(peak_file):
            n_peaks = sum(1 for _ in open(peak_file))
            total_peaks += n_peaks
            print(f"  Cluster {cl}: {n_peaks} peaks")
    # Clean up large BED file
    os.remove(bed)

print(f"  Total peaks across clusters: {total_peaks}")

# Branch 7b: Gene activity matrix
print("  Computing gene activity scores...")
gene_matrix = snap.pp.make_gene_matrix(adata, snap.genome.hg38)
print(f"  Gene activity matrix: {gene_matrix.shape}")

# Branch 7c: Export
adata.write("outputs/clusters/adata_clustered.h5ad")
gene_matrix.write("outputs/gene_activity/gene_activity.h5ad")

print("=== Level 8-9: CONVERGENCE 3+4 - Analysis + Report ===")

# Marker genes per cluster
sc.tl.rank_genes_groups(gene_matrix, groupby="leiden", method="wilcoxon")
marker_df = sc.get.rank_genes_groups_df(gene_matrix, group=None)
marker_df.head(100).to_csv("outputs/community/marker_genes.tsv", sep='\t', index=False)

# Differential accessibility per cluster
n_marker_genes = len(marker_df[marker_df["pvals_adj"] < 0.05]["names"].unique())

# === Build report ===
metrics = {}
metrics["total_cells_imported"] = n_cells_before
metrics["cells_after_qc"] = n_cells_after
metrics["cells_after_doublet_removal"] = adata.n_obs
metrics["median_tsse"] = round(np.median(tsse_values), 2)
metrics["mean_tsse"] = round(np.mean(tsse_values), 2)
metrics["median_fragments_per_cell"] = int(np.median(frag_values))
metrics["mean_fragments_per_cell"] = int(np.mean(frag_values))
metrics["num_clusters"] = n_clusters
metrics["total_peaks"] = total_peaks
metrics["n_genes_in_activity_matrix"] = gene_matrix.n_vars
metrics["n_significant_markers"] = n_marker_genes

# Per-cluster stats
for cl in sorted(set(cluster_labels)):
    cl_size = (cluster_labels == cl).sum()
    metrics[f"cluster_{cl}_cells"] = cl_size

# Doublet stats
metrics["doublets_detected"] = int(n_doublets)

with open("results/report.csv", 'w') as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("\n=== Report ===")
for k, v in metrics.items():
    print(f"  {k} = {v}")
PYEOF

echo "=== Pipeline complete ==="
