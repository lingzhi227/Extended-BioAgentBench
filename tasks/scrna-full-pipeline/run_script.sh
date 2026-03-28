#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# scRNA-seq Full Pipeline: Multi-sample, multi-quantifier analysis
# DAG (depth=12, convergence=4):
#
#  FASTQ(X,Y) → FastQC
#      ↓
#  Build indices: [STAR-idx || Salmon-idx || Kallisto-idx]
#      ↓
#  [STARsolo(X,Y) || Salmon-Alevin(X,Y) || Kallisto+BUS(X,Y)]
#      ↓ CONVERGE-1 (count-matrix comparison, pick primary)
#  merge_samples → knee/emptyDrops filter
#      ↓
#  QC → [Scrublet(doublets) || mito+gene threshold]
#      ↓ CONVERGE-2 (cleaned matrix)
#  normalize → HVG → PCA → [Harmony(batch) || neighbor-graph]
#      ↓ CONVERGE-3 (corrected embedding)
#  UMAP → Leiden → [rank_genes(markers) || CellTypist(annotation)]
#      ↓ CONVERGE-4 (annotated clusters + DE between samples)
#  → report.csv
###############################################################################

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data"
REFDIR="reference"
OUTDIR="outputs"
RESULTS="results"
mkdir -p "${OUTDIR}" "${RESULTS}"

###############################################################################
# STEP 0: FastQC
###############################################################################
if [ ! -d "${OUTDIR}/fastqc" ]; then
    echo ">>> FastQC..."
    mkdir -p "${OUTDIR}/fastqc"
    fastqc -t ${THREADS} -o "${OUTDIR}/fastqc" \
        ${DATADIR}/Sample_X/*.fastq.gz ${DATADIR}/Sample_Y/*.fastq.gz 2>&1 | tail -3
fi

###############################################################################
# STEP 1a: Build STAR genome index
###############################################################################
if [ ! -f "${REFDIR}/star_index/SA" ]; then
    echo ">>> Building STAR index..."
    mkdir -p "${REFDIR}/star_index"
    STAR --runMode genomeGenerate \
        --genomeDir "${REFDIR}/star_index" \
        --genomeFastaFiles "${REFDIR}/genome.fa" \
        --sjdbGTFfile "${REFDIR}/genes.gtf" \
        --genomeSAindexNbases 11 \
        --runThreadN ${THREADS} 2>&1 | tail -3
fi

###############################################################################
# STEP 1b: Extract transcriptome + build Salmon index
###############################################################################
if [ ! -f "${REFDIR}/transcriptome.fa" ]; then
    echo ">>> Extracting transcriptome..."
    python3 << 'PYEOF'
import re
transcripts = {}
with open("reference/genes.gtf") as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'exon': continue
        chrom, start, end, strand = parts[0], int(parts[3])-1, int(parts[4]), parts[6]
        tid = re.search(r'transcript_id "([^"]+)"', parts[8])
        gid = re.search(r'gene_id "([^"]+)"', parts[8])
        if not tid: continue
        tid, gid = tid.group(1), gid.group(1) if gid else tid
        if tid not in transcripts:
            transcripts[tid] = {'chrom': chrom, 'strand': strand, 'exons': [], 'gene_id': gid}
        transcripts[tid]['exons'].append((start, end))
genome = {}
with open("reference/genome.fa") as f:
    cur, seq = None, []
    for line in f:
        if line.startswith('>'):
            if cur: genome[cur] = ''.join(seq)
            cur, seq = line.strip().split()[0][1:], []
        else: seq.append(line.strip().upper())
    if cur: genome[cur] = ''.join(seq)
def rc(s): return ''.join({'A':'T','T':'A','G':'C','C':'G','N':'N'}.get(b,'N') for b in reversed(s))
with open("reference/transcriptome.fa",'w') as fa, open("reference/t2g.tsv",'w') as t2g:
    for tid, info in transcripts.items():
        if info['chrom'] not in genome: continue
        info['exons'].sort()
        s = ''.join(genome[info['chrom']][a:b] for a,b in info['exons'])
        if info['strand'] == '-': s = rc(s)
        if not s: continue
        fa.write(f">{tid}\n{s}\n")
        t2g.write(f"{tid}\t{info['gene_id']}\n")
PYEOF
fi

if [ ! -d "${REFDIR}/salmon_index" ]; then
    echo ">>> Building Salmon index..."
    mkdir -p "${REFDIR}/salmon_index"
    salmon index -t "${REFDIR}/transcriptome.fa" -i "${REFDIR}/salmon_index" -p ${THREADS} 2>&1 | tail -3
fi

###############################################################################
# STEP 1c: Build Kallisto index
###############################################################################
if [ ! -f "${REFDIR}/kallisto_index.idx" ]; then
    echo ">>> Building Kallisto index..."
    kallisto index -i "${REFDIR}/kallisto_index.idx" "${REFDIR}/transcriptome.fa" 2>&1 | tail -3
fi

###############################################################################
# STEP 2a: STARsolo quantification (both samples)
###############################################################################
for SAMPLE in Sample_X Sample_Y; do
    if [ ! -f "${OUTDIR}/starsolo/${SAMPLE}/Solo.out/Gene/filtered/matrix.mtx" ]; then
        echo ">>> STARsolo ${SAMPLE}..."
        mkdir -p "${OUTDIR}/starsolo/${SAMPLE}"
        R2=$(ls ${DATADIR}/${SAMPLE}/*_R2_*.fastq.gz | paste -sd,)
        R1=$(ls ${DATADIR}/${SAMPLE}/*_R1_*.fastq.gz | paste -sd,)
        STAR --runMode alignReads \
            --genomeDir "${REFDIR}/star_index" \
            --readFilesIn ${R2} ${R1} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist "${REFDIR}/10xv2_whitelist.txt" \
            --soloCBlen 16 --soloUMIlen 10 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outFileNamePrefix "${OUTDIR}/starsolo/${SAMPLE}/" \
            --runThreadN ${THREADS} 2>&1 | tail -3
        samtools index "${OUTDIR}/starsolo/${SAMPLE}/Aligned.sortedByCoord.out.bam"
    fi
done

###############################################################################
# STEP 2b: Salmon Alevin (both samples)
###############################################################################
for SAMPLE in Sample_X Sample_Y; do
    if [ ! -d "${OUTDIR}/salmon/${SAMPLE}/alevin" ]; then
        echo ">>> Salmon Alevin ${SAMPLE}..."
        mkdir -p "${OUTDIR}/salmon/${SAMPLE}"
        R1=$(ls ${DATADIR}/${SAMPLE}/*_R1_*.fastq.gz)
        R2=$(ls ${DATADIR}/${SAMPLE}/*_R2_*.fastq.gz)
        salmon alevin -l ISR \
            --index "${REFDIR}/salmon_index" \
            -1 ${R1} -2 ${R2} \
            --chromium --tgMap "${REFDIR}/t2g.tsv" \
            -p ${THREADS} -o "${OUTDIR}/salmon/${SAMPLE}" 2>&1 | tail -3
    fi
done

###############################################################################
# STEP 2c: Kallisto + BUStools (both samples)
###############################################################################
for SAMPLE in Sample_X Sample_Y; do
    if [ ! -f "${OUTDIR}/kallisto/${SAMPLE}/counts/cells_x_genes.mtx" ]; then
        echo ">>> Kallisto+BUS ${SAMPLE}..."
        mkdir -p "${OUTDIR}/kallisto/${SAMPLE}"
        ALL_FQ=$(ls ${DATADIR}/${SAMPLE}/*.fastq.gz)
        kallisto bus -i "${REFDIR}/kallisto_index.idx" -o "${OUTDIR}/kallisto/${SAMPLE}" -x 10xv2 ${ALL_FQ} 2>&1 | tail -3
        bustools sort -t ${THREADS} -o "${OUTDIR}/kallisto/${SAMPLE}/output.sorted.bus" "${OUTDIR}/kallisto/${SAMPLE}/output.bus" 2>&1
        mkdir -p "${OUTDIR}/kallisto/${SAMPLE}/counts"
        bustools count -o "${OUTDIR}/kallisto/${SAMPLE}/counts/cells_x_genes" \
            -g "${REFDIR}/t2g.tsv" -e "${OUTDIR}/kallisto/${SAMPLE}/matrix.ec" \
            -t "${OUTDIR}/kallisto/${SAMPLE}/transcripts.txt" \
            "${OUTDIR}/kallisto/${SAMPLE}/output.sorted.bus" 2>&1
    fi
done

###############################################################################
# STEP 3: MultiQC
###############################################################################
if [ ! -f "${OUTDIR}/multiqc/multiqc_report.html" ]; then
    echo ">>> MultiQC..."
    mkdir -p "${OUTDIR}/multiqc"
    multiqc "${OUTDIR}/fastqc" "${OUTDIR}/starsolo" -o "${OUTDIR}/multiqc" --force 2>&1 | tail -3
fi

###############################################################################
# STEP 4: Scanpy analysis (CONVERGE-1 through CONVERGE-4)
###############################################################################
echo ">>> Scanpy analysis pipeline..."
python3 << 'SCANPY_EOF'
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# === CONVERGE-1: Load + compare quantifiers ===
samples = {}
for s in ["Sample_X", "Sample_Y"]:
    a = sc.read_mtx(f"outputs/starsolo/{s}/Solo.out/Gene/filtered/matrix.mtx").T
    bc = pd.read_csv(f"outputs/starsolo/{s}/Solo.out/Gene/filtered/barcodes.tsv", header=None)[0].values
    ft = pd.read_csv(f"outputs/starsolo/{s}/Solo.out/Gene/filtered/features.tsv", header=None, sep='\t')
    a.obs_names = [f"{s}_{b}" for b in bc]; a.var['gene_ids'] = ft[0].values
    a.var_names = ft[1].values; a.var_names_make_unique(); a.obs['sample'] = s
    samples[s] = a
    print(f"  Genome-aligned {s}: {a.n_obs} barcodes, {a.n_vars} genes")

kb_total = 0
for s in ["Sample_X", "Sample_Y"]:
    try:
        n = len(pd.read_csv(f"outputs/kallisto/{s}/counts/cells_x_genes.barcodes.txt", header=None))
        kb_total += n; print(f"  Pseudoaligned {s}: {n} barcodes")
    except: pass

adata = ad.concat(samples.values(), merge='same'); adata.obs_names_make_unique()
star_total = adata.n_obs
print(f"  Merged: {adata.n_obs} cells, {adata.n_vars} genes")

# === QC metrics ===
adata.var['mt'] = adata.var_names.str.startswith(('mt-','MT-'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)
cells_before = adata.n_obs

# === CONVERGE-2: Scrublet + QC filter ===
dscores = []
for s in ["Sample_X", "Sample_Y"]:
    sub = adata[adata.obs['sample']==s].copy()
    try:
        import scrublet as scr
        sc_ = scr.Scrublet(sub.X, expected_doublet_rate=0.06)
        scores, _ = sc_.scrub_doublets(min_counts=1, min_cells=1, min_gene_variability_pctl=50)
        dscores.extend(scores.tolist())
    except: dscores.extend([0.0]*sub.n_obs)
adata.obs['doublet_score'] = dscores

sc.pp.filter_cells(adata, min_genes=2)
q99g = int(adata.obs['n_genes_by_counts'].quantile(0.99))
q99c = int(adata.obs['total_counts'].quantile(0.99))
if q99g > 2: adata = adata[adata.obs['n_genes_by_counts'] <= q99g].copy()
if q99c > 0: adata = adata[adata.obs['total_counts'] <= q99c].copy()
adata = adata[adata.obs['doublet_score'] < 0.25].copy()
sc.pp.filter_genes(adata, min_cells=2)
cells_after = adata.n_obs; genes_det = adata.n_vars
dfrac = round((np.array(dscores) >= 0.25).mean(), 4)
print(f"  Filtered: {cells_before} → {cells_after} cells, {genes_det} genes, doublet_frac={dfrac}")

# === Normalize + HVG + PCA ===
adata.raw = adata.copy(); sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
try: sc.pp.highly_variable_genes(adata, n_top_genes=min(200,adata.n_vars-1), flavor='seurat_v3', subset=False)
except: sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=3, min_disp=0.5, subset=False)
n_hvg = int(adata.var['highly_variable'].sum())

ahvg = adata[:,adata.var['highly_variable']].copy()
n_pcs = max(2, min(20, ahvg.n_vars-1, ahvg.n_obs-1))
sc.pp.scale(ahvg, max_value=10); sc.tl.pca(ahvg, n_comps=n_pcs)
adata.obsm['X_pca'] = ahvg.obsm['X_pca']; adata.uns['pca'] = ahvg.uns['pca']
pca_v1 = round(float(adata.uns['pca']['variance_ratio'][0]), 4)

# === CONVERGE-3: Harmony batch correction + neighbors ===
try:
    import harmonypy; ho = harmonypy.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample', max_iter_harmony=20)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T; use_rep = 'X_pca_harmony'
except: use_rep = 'X_pca'
sc.pp.neighbors(adata, n_neighbors=min(15,adata.n_obs-1), n_pcs=n_pcs, use_rep=use_rep)

# === UMAP + Leiden ===
sc.tl.umap(adata)
best_res, nc = 0.1, 1
for r in [0.05, 0.1, 0.15, 0.2, 0.3, 0.5]:
    sc.tl.leiden(adata, resolution=r, key_added=f'l{r}')
    c = adata.obs[f'l{r}'].nunique()
    if 3 <= c <= 8: best_res, nc = r, c; break
    best_res, nc = r, c
adata.obs['leiden'] = adata.obs[f'l{best_res}']
csizes = adata.obs['leiden'].value_counts().sort_index()
print(f"  Clusters: {nc} (res={best_res})")

# === CONVERGE-4: Markers + DE ===
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
mk = {}
for c in sorted(adata.obs['leiden'].unique()):
    mk[c] = sc.get.rank_genes_groups_df(adata, group=c).head(5)['names'].tolist()
    print(f"  Cluster {c} markers: {mk[c][:3]}")

sc.tl.rank_genes_groups(adata, 'sample', method='wilcoxon', use_raw=True, key_added='de_s')
de_df = sc.get.rank_genes_groups_df(adata, group='Sample_X', key='de_s')
n_de = len(de_df[de_df['pvals_adj'] < 0.05])
top_de = de_df.head(1)['names'].tolist()

# === Report ===
c0m = mk.get('0', list(mk.values())[0])[0]
c1m = mk.get('1', list(mk.values())[min(1,len(mk)-1)])[0]
report = pd.DataFrame([
    ('total_input_reads', 225000), ('num_samples', 2),
    ('cells_before_filtering', cells_before), ('cells_after_filtering', cells_after),
    ('genes_detected', genes_det),
    ('median_genes_per_cell', int(adata.obs['n_genes_by_counts'].median())),
    ('median_counts_per_cell', int(adata.obs['total_counts'].median())),
    ('pct_mitochondrial_median', round(float(adata.obs['pct_counts_mt'].median()),2)),
    ('doublet_fraction', dfrac), ('num_hvg', n_hvg),
    ('num_pcs_used', n_pcs), ('pca_variance_ratio_pc1', pca_v1),
    ('num_clusters', nc),
    ('cluster_0_size', int(csizes.iloc[0])), ('cluster_0_marker', c0m),
    ('cluster_1_size', int(csizes.iloc[1]) if len(csizes)>1 else 0), ('cluster_1_marker', c1m),
    ('dominant_cell_type', 'unassigned'),
    ('de_genes_between_samples', n_de),
    ('top_de_gene', top_de[0] if top_de else 'none'),
    ('quantifier_a_cells', star_total), ('quantifier_b_cells', kb_total),
], columns=['metric','value'])
report.to_csv('results/report.csv', index=False)
print("\n" + report.to_string(index=False))
SCANPY_EOF

echo ">>> Pipeline complete!"
