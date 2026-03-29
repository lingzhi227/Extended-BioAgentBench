#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Spatial Transcriptomics (Visium) Analysis Pipeline
# DAG (depth=11, convergence=4):
#
#  [FASTQ → FastQC] || [BAM → samtools flagstat+idxstats]
#      ↓ CONVERGE-1 (read QC + alignment QC)
#  Load Visium data (matrix + spatial + image) → filter → normalize
#      ↓
#  [squidpy spatial-autocorrelation || Leiden clustering]
#      ↓ CONVERGE-2 (spatial vs non-spatial clustering)
#  [squidpy neighborhood-enrichment || squidpy co-occurrence]
#      ↓ CONVERGE-3 (neighborhood metrics)
#  [rank_genes_groups(markers) || squidpy spatial-variable]
#      ↓ CONVERGE-4 (markers + spatial genes)
#  MultiQC → report.csv
###############################################################################

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data"
OUTDIR="outputs"
RESULTS="results"
mkdir -p "${OUTDIR}" "${RESULTS}"

###############################################################################
# STEP 0: FastQC on raw Visium reads
###############################################################################
if [ ! -d "${OUTDIR}/fastqc" ]; then
    echo ">>> FastQC..."
    mkdir -p "${OUTDIR}/fastqc"
    fastqc -t ${THREADS} -o "${OUTDIR}/fastqc" \
        ${DATADIR}/R1_L001.fastq.gz ${DATADIR}/R2_L001.fastq.gz \
        ${DATADIR}/R1_L002.fastq.gz ${DATADIR}/R2_L002.fastq.gz \
        2>&1 | tail -3
fi

###############################################################################
# STEP 1: Alignment QC from pre-aligned BAM
###############################################################################
if [ ! -f "${OUTDIR}/bam_qc/flagstat.txt" ]; then
    echo ">>> samtools QC..."
    mkdir -p "${OUTDIR}/bam_qc"
    samtools flagstat "${DATADIR}/outs/possorted_genome_bam.bam" > "${OUTDIR}/bam_qc/flagstat.txt" 2>/dev/null
    samtools idxstats "${DATADIR}/outs/possorted_genome_bam.bam" > "${OUTDIR}/bam_qc/idxstats.txt" 2>/dev/null
    samtools stats "${DATADIR}/outs/possorted_genome_bam.bam" > "${OUTDIR}/bam_qc/stats.txt" 2>/dev/null || true
fi

###############################################################################
# STEP 2: MultiQC aggregation
###############################################################################
if [ ! -f "${OUTDIR}/multiqc/multiqc_report.html" ]; then
    echo ">>> MultiQC..."
    mkdir -p "${OUTDIR}/multiqc"
    multiqc "${OUTDIR}/fastqc" "${OUTDIR}/bam_qc" \
        -o "${OUTDIR}/multiqc" --force 2>&1 | tail -3
fi

###############################################################################
# STEP 3: Spatial analysis pipeline
###############################################################################
echo ">>> Spatial analysis..."
python3 << 'SPATIAL_EOF'
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# === CONVERGE-1: Load alignment QC + Visium data ===
print("=== Loading data ===")
with open("outputs/bam_qc/flagstat.txt") as f:
    flagstat = f.read()
total_reads = int(flagstat.split('\n')[0].split('+')[0].strip())
mapped_reads = int([l for l in flagstat.split('\n') if 'primary mapped' in l][0].split('+')[0].strip())
dup_reads = int([l for l in flagstat.split('\n') if 'duplicates' in l][0].split('+')[0].strip())
mapping_rate = round(mapped_reads/total_reads*100, 2) if total_reads > 0 else 0
print(f"  BAM: {total_reads} reads, {mapping_rate}% mapped, {dup_reads} dups")

adata = sc.read_visium("data/outs/", count_file="raw_feature_bc_matrix.h5")
# Use gene symbols instead of Ensembl IDs
if 'gene_ids' not in adata.var.columns and adata.var_names[0].startswith('ENSG'):
    adata.var['gene_ids'] = adata.var_names.copy()
    # Read features.tsv for gene names
    import gzip
    feat = pd.read_csv("data/outs/raw_feature_bc_matrix/features.tsv.gz", sep='\t', header=None)
    id_to_name = dict(zip(feat[0], feat[1]))
    adata.var_names = [id_to_name.get(g, g) for g in adata.var_names]
    adata.var_names_make_unique()
print(f"  Visium: {adata.n_obs} spots, {adata.n_vars} genes, {adata.X.sum():.0f} counts")

# Spatial extent
spatial = adata.obsm['spatial']
x_range = round(float(spatial[:,0].max() - spatial[:,0].min()), 1)
y_range = round(float(spatial[:,1].max() - spatial[:,1].min()), 1)

# In-tissue spots
pos = pd.read_csv("data/outs/spatial/tissue_positions.csv")
in_tissue = int((pos.iloc[:,1]==1).sum())
print(f"  In-tissue spots: {in_tissue}")

# === QC + filter ===
print("\n=== QC ===")
sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)
spots_before = adata.n_obs
sc.pp.filter_cells(adata, min_counts=1)
sc.pp.filter_genes(adata, min_cells=1)
spots_after = adata.n_obs
genes_det = adata.n_vars
med_genes = int(adata.obs['n_genes_by_counts'].median()) if spots_after > 0 else 0
med_counts = int(adata.obs['total_counts'].median()) if spots_after > 0 else 0
print(f"  Filtered: {spots_before} → {spots_after} spots, {genes_det} genes")

# === Normalize + HVG + PCA ===
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
n_top = min(100, adata.n_vars-1)
try:
    sc.pp.highly_variable_genes(adata, n_top_genes=max(2,n_top), flavor='seurat_v3', subset=False)
except:
    adata.var['highly_variable'] = True
n_hvg = int(adata.var['highly_variable'].sum())
ahvg = adata[:,adata.var['highly_variable']].copy()
n_pcs = max(2, min(15, ahvg.n_vars-1, ahvg.n_obs-1))
sc.pp.scale(ahvg, max_value=10); sc.tl.pca(ahvg, n_comps=n_pcs)
adata.obsm['X_pca'] = ahvg.obsm['X_pca']; adata.uns['pca'] = ahvg.uns['pca']
pca_v1 = round(float(adata.uns['pca']['variance_ratio'][0]), 4) if spots_after > 2 else 0

# === CONVERGE-2: Spatial autocorrelation + Leiden ===
print("\n=== Clustering ===")
nn = min(10, adata.n_obs-1)
sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=n_pcs); sc.tl.umap(adata)
best_res, nc = 0.3, 1
for r in [0.1, 0.2, 0.3, 0.5, 0.8]:
    sc.tl.leiden(adata, resolution=r, key_added=f'l{r}')
    c = adata.obs[f'l{r}'].nunique()
    if 2 <= c <= 8: best_res, nc = r, c; break
    best_res, nc = r, c
adata.obs['leiden'] = adata.obs[f'l{best_res}']
print(f"  Leiden: {nc} clusters (res={best_res})")

morans_genes, top_mi, n_sig = [], 0, 0
try:
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
    sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=1)
    mi = adata.uns['moranI'].sort_values('I', ascending=False)
    morans_genes = mi.head(5).index.tolist()
    top_mi = round(float(mi['I'].iloc[0]), 4)
    n_sig = len(mi[mi['pval_norm'] < 0.05])
    print(f"  Moran's I: top={morans_genes[0]} (I={top_mi}), {n_sig} significant")
except Exception as e:
    print(f"  Moran's I: {e}")

# === CONVERGE-3: Neighborhood enrichment + Co-occurrence ===
print("\n=== Neighborhood ===")
max_nz, max_co = 0, 0
try:
    sq.gr.nhood_enrichment(adata, cluster_key='leiden')
    max_nz = round(float(np.nanmax(adata.uns['leiden_nhood_enrichment']['zscore'])), 4)
    print(f"  Nhood max z: {max_nz}")
except Exception as e:
    print(f"  Nhood: {e}")
try:
    sq.gr.co_occurrence(adata, cluster_key='leiden')
    max_co = round(float(np.nanmax(adata.uns['leiden_co_occurrence']['occ'])), 4)
    print(f"  Co-occ max: {max_co}")
except Exception as e:
    print(f"  Co-occ: {e}")

# === CONVERGE-4: Markers + spatial variable ===
print("\n=== Markers ===")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
mk = {}
for c in sorted(adata.obs['leiden'].unique()):
    mk[c] = sc.get.rank_genes_groups_df(adata, group=c).head(3)['names'].tolist()
    print(f"  Cluster {c}: {mk[c][:2]}")

# === Report ===
csizes = adata.obs['leiden'].value_counts().sort_index()
c0m = mk.get('0', list(mk.values())[0])[0]
report = pd.DataFrame([
    ('total_reads', total_reads), ('mapped_reads', mapped_reads),
    ('mapping_rate_pct', mapping_rate), ('duplicate_reads', dup_reads),
    ('total_spots', spots_before), ('in_tissue_spots', in_tissue),
    ('spots_with_counts', spots_after), ('genes_detected', genes_det),
    ('median_genes_per_spot', med_genes), ('median_counts_per_spot', med_counts),
    ('num_hvg', n_hvg), ('num_pcs', n_pcs), ('pca_variance_ratio_pc1', pca_v1),
    ('num_spatial_clusters', nc),
    ('cluster_0_size', int(csizes.iloc[0]) if len(csizes) > 0 else 0),
    ('cluster_0_marker', c0m),
    ('top_spatial_gene', morans_genes[0] if morans_genes else 'none'),
    ('top_morans_I', top_mi), ('spatially_variable_genes', n_sig),
    ('nhood_enrichment_max_z', max_nz), ('co_occurrence_max', max_co),
    ('spatial_extent_x', x_range), ('spatial_extent_y', y_range),
], columns=['metric','value'])
report.to_csv('results/report.csv', index=False)
print("\n" + report.to_string(index=False))
SPATIAL_EOF

echo ">>> Pipeline complete!"
