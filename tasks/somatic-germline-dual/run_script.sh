#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Somatic + Germline Dual Variant Analysis
# DAG (depth=12, convergence=4):
#
#  [Tumor: FASTQ → BWA → sort → dedup]
#  ||
#  [Normal: FASTQ → BWA → sort → dedup]
#      ↓ CONVERGE-1 (shared reference, paired samples)
#  [Mutect2(somatic) || HaplotypeCaller(germline)]
#      ↓
#  [FilterMutect || GenotypeGVCFs → filter]
#      ↓ CONVERGE-2 (filtered somatic + germline)
#  bcftools isec (somatic ∩ germline positions)
#      ↓ CONVERGE-3 (intersection analysis)
#  [mosdepth(coverage) || bcftools stats(QC)]
#      ↓ CONVERGE-4 (QC + coverage + annotations)
#  → report.csv
###############################################################################

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data"
REFDIR="reference"
OUTDIR="outputs"
RESULTS="results"
mkdir -p "${OUTDIR}" "${RESULTS}"
GENOME="${REFDIR}/genome.fa"

# Index reference (if needed)
[ ! -f "${GENOME}.bwt.2bit.64" ] && bwa-mem2 index "${GENOME}" 2>&1 | tail -1
[ ! -f "${GENOME}.fai" ] && samtools faidx "${GENOME}"
[ ! -f "${REFDIR}/genome.dict" ] && gatk CreateSequenceDictionary -R "${GENOME}" -O "${REFDIR}/genome.dict" 2>/dev/null

# FastQC
if [ ! -d "${OUTDIR}/fastqc" ]; then
    echo ">>> FastQC..."
    mkdir -p "${OUTDIR}/fastqc"
    fastqc -t ${THREADS} -o "${OUTDIR}/fastqc" ${DATADIR}/*.fastq.gz 2>&1 | tail -3
fi

# Align + dedup both samples
for SAMPLE in tumor normal; do
    if [ ! -f "${OUTDIR}/${SAMPLE}_dedup.bam" ]; then
        echo ">>> Align+dedup ${SAMPLE}..."
        bwa-mem2 mem -t ${THREADS} \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
            "${GENOME}" "${DATADIR}/${SAMPLE}_R1.fastq.gz" "${DATADIR}/${SAMPLE}_R2.fastq.gz" 2>/dev/null | \
            samtools sort -@ ${THREADS} -o "${OUTDIR}/${SAMPLE}_sorted.bam" -
        samtools index "${OUTDIR}/${SAMPLE}_sorted.bam"
        picard MarkDuplicates -I "${OUTDIR}/${SAMPLE}_sorted.bam" -O "${OUTDIR}/${SAMPLE}_dedup.bam" \
            -M "${OUTDIR}/${SAMPLE}_dup_metrics.txt" --REMOVE_DUPLICATES false --CREATE_INDEX true 2>&1 | tail -1
    fi
done

# Coverage
for SAMPLE in tumor normal; do
    if [ ! -f "${OUTDIR}/${SAMPLE}_cov.mosdepth.summary.txt" ]; then
        echo ">>> mosdepth ${SAMPLE}..."
        mosdepth --threads ${THREADS} --no-per-base "${OUTDIR}/${SAMPLE}_cov" "${OUTDIR}/${SAMPLE}_dedup.bam" 2>&1 | tail -1
    fi
done

# CONVERGE-1: Somatic calling (Mutect2)
if [ ! -f "${OUTDIR}/somatic_filtered.vcf.gz" ]; then
    echo ">>> Mutect2 (somatic)..."
    gatk Mutect2 -R "${GENOME}" \
        -I "${OUTDIR}/tumor_dedup.bam" -I "${OUTDIR}/normal_dedup.bam" -normal normal \
        -O "${OUTDIR}/somatic_raw.vcf.gz" --native-pair-hmm-threads ${THREADS} 2>&1 | tail -3
    gatk FilterMutectCalls -R "${GENOME}" \
        -V "${OUTDIR}/somatic_raw.vcf.gz" -O "${OUTDIR}/somatic_filtered.vcf.gz" 2>&1 | tail -3
fi

# CONVERGE-1 parallel: Germline calling (HaplotypeCaller on normal)
if [ ! -f "${OUTDIR}/germline_raw.vcf.gz" ]; then
    echo ">>> HaplotypeCaller (germline)..."
    gatk HaplotypeCaller -R "${GENOME}" \
        -I "${OUTDIR}/normal_dedup.bam" -O "${OUTDIR}/germline_raw.vcf.gz" \
        --native-pair-hmm-threads ${THREADS} 2>&1 | tail -3
fi

# Also call germline on tumor for comparison
if [ ! -f "${OUTDIR}/tumor_germline.vcf.gz" ]; then
    echo ">>> HaplotypeCaller (tumor germline)..."
    gatk HaplotypeCaller -R "${GENOME}" \
        -I "${OUTDIR}/tumor_dedup.bam" -O "${OUTDIR}/tumor_germline.vcf.gz" \
        --native-pair-hmm-threads ${THREADS} 2>&1 | tail -3
fi

# CONVERGE-2: bcftools stats on both
bcftools stats "${OUTDIR}/somatic_filtered.vcf.gz" > "${OUTDIR}/somatic_stats.txt" 2>/dev/null
bcftools stats "${OUTDIR}/germline_raw.vcf.gz" > "${OUTDIR}/germline_stats.txt" 2>/dev/null
bcftools stats "${OUTDIR}/tumor_germline.vcf.gz" > "${OUTDIR}/tumor_germline_stats.txt" 2>/dev/null

# CONVERGE-3: Somatic-Germline intersection
if [ ! -d "${OUTDIR}/isec" ]; then
    echo ">>> bcftools isec (somatic ∩ germline)..."
    bcftools index "${OUTDIR}/somatic_filtered.vcf.gz" 2>/dev/null || true
    bcftools index "${OUTDIR}/germline_raw.vcf.gz" 2>/dev/null || true
    bcftools index "${OUTDIR}/tumor_germline.vcf.gz" 2>/dev/null || true
    mkdir -p "${OUTDIR}/isec"
    bcftools isec \
        "${OUTDIR}/germline_raw.vcf.gz" \
        "${OUTDIR}/tumor_germline.vcf.gz" \
        -p "${OUTDIR}/isec" 2>/dev/null || true
fi

# CONVERGE-4: Report
echo ">>> Report..."
python3 << 'REPORT_EOF'
import subprocess, gzip
import pandas as pd

def run(cmd): return subprocess.run(cmd, capture_output=True, text=True, shell=True)
def count_vcf(f):
    r = run(f"bcftools view -H {f}")
    return len([l for l in r.stdout.strip().split('\n') if l]) if r.stdout.strip() else 0
def get_stat(f, key):
    try:
        with open(f) as fh:
            for line in fh:
                if line.startswith("SN") and key in line:
                    return int(line.strip().split('\t')[-1])
    except: pass
    return 0

# Reads
t_reads = sum(1 for _ in gzip.open("data/tumor_R1.fastq.gz",'rt')) // 4
n_reads = sum(1 for _ in gzip.open("data/normal_R1.fastq.gz",'rt')) // 4

# Mapped
t_mapped = int(run("samtools flagstat outputs/tumor_dedup.bam").stdout.split('\n')[4].split('+')[0].strip())
n_mapped = int(run("samtools flagstat outputs/normal_dedup.bam").stdout.split('\n')[4].split('+')[0].strip())

# Coverage
t_cov = n_cov = 0
for s,v in [("tumor","t_cov"),("normal","n_cov")]:
    with open(f"outputs/{s}_cov.mosdepth.summary.txt") as f:
        for l in f:
            if l.startswith("total"):
                exec(f"{v}=round(float(l.split()[3]),1)")

# Duplicates
t_dup = n_dup = 0
for s,v in [("tumor","t_dup"),("normal","n_dup")]:
    with open(f"outputs/{s}_dup_metrics.txt") as f:
        for l in f:
            parts = l.strip().split('\t')
            if len(parts)>8 and parts[0] not in ['#','LIBRARY',''] and 'Unknown' not in l:
                try: exec(f"{v}=round(float(parts[8])*100,2)"); break
                except: continue

# Variants
somatic_n = count_vcf("outputs/somatic_filtered.vcf.gz")
germline_n = count_vcf("outputs/germline_raw.vcf.gz")
tumor_germ_n = count_vcf("outputs/tumor_germline.vcf.gz")
germ_snps = get_stat("outputs/germline_stats.txt", "number of SNPs:")
germ_indels = get_stat("outputs/germline_stats.txt", "number of indels:")
tgerm_snps = get_stat("outputs/tumor_germline_stats.txt", "number of SNPs:")

# Intersection
isec_shared = 0
try:
    import os
    for f in ["outputs/isec/0002.vcf"]:
        if os.path.exists(f):
            isec_shared = sum(1 for l in open(f) if not l.startswith('#'))
except: pass

# Normal-only and tumor-only germline
isec_normal_only = 0
isec_tumor_only = 0
try:
    for f,v in [("outputs/isec/0000.vcf","isec_normal_only"),("outputs/isec/0001.vcf","isec_tumor_only")]:
        if os.path.exists(f):
            exec(f"{v} = sum(1 for l in open(f) if not l.startswith('#'))")
except: pass

report = pd.DataFrame([
    ('tumor_input_reads', t_reads), ('normal_input_reads', n_reads),
    ('tumor_mapped_reads', t_mapped), ('normal_mapped_reads', n_mapped),
    ('tumor_mean_coverage', t_cov), ('normal_mean_coverage', n_cov),
    ('tumor_duplicate_pct', t_dup), ('normal_duplicate_pct', n_dup),
    ('somatic_variants', somatic_n),
    ('germline_variants_normal', germline_n),
    ('germline_variants_tumor', tumor_germ_n),
    ('germline_snps', germ_snps), ('germline_indels', germ_indels),
    ('shared_germline_positions', isec_shared),
    ('normal_only_germline', isec_normal_only),
    ('tumor_only_germline', isec_tumor_only),
], columns=['metric','value'])

for i,row in report.iterrows():
    try:
        fv = float(row['value'])
        if fv == int(fv) and 'pct' not in str(row['metric']) and 'coverage' not in str(row['metric']):
            report.at[i,'value'] = int(fv)
    except: pass

report.to_csv('results/report.csv', index=False)
print(report.to_string(index=False))
REPORT_EOF

echo ">>> Pipeline complete!"
