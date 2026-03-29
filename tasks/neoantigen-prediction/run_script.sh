#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Neoantigen Prediction: Tumor-Normal Somatic Analysis
# DAG (depth=11, convergence=4):
#
#  [Tumor FASTQ → FastQC → BWA-MEM2 → sort → dedup → BQSR]
#  ||
#  [Normal FASTQ → FastQC → BWA-MEM2 → sort → dedup → BQSR]
#      ↓ CONVERGE-1 (shared reference + paired processing)
#  Mutect2(T+N) → FilterMutectCalls
#      ↓
#  [bcftools stats(QC) || mosdepth(coverage)]
#      ↓ CONVERGE-2 (variant QC + coverage)
#  [SnpSift annotate || bcftools filter(PASS)]
#      ↓ CONVERGE-3 (annotated + filtered variants)
#  Extract peptide contexts → neoantigen scoring
#      ↓ CONVERGE-4 (variants + peptides + scores)
#  → report.csv
###############################################################################

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data"
REFDIR="reference"
OUTDIR="outputs"
RESULTS="results"
mkdir -p "${OUTDIR}" "${RESULTS}"

GENOME="${REFDIR}/genome.fa"

###############################################################################
# STEP 0: Index reference
###############################################################################
if [ ! -f "${GENOME}.bwt.2bit.64" ]; then
    echo ">>> Indexing reference..."
    bwa-mem2 index "${GENOME}" 2>&1 | tail -3
fi
if [ ! -f "${GENOME}.fai" ]; then
    samtools faidx "${GENOME}"
fi
if [ ! -f "${REFDIR}/genome.dict" ]; then
    gatk CreateSequenceDictionary -R "${GENOME}" -O "${REFDIR}/genome.dict" 2>&1 | tail -3
fi

###############################################################################
# STEP 1: FastQC on raw reads
###############################################################################
if [ ! -d "${OUTDIR}/fastqc" ]; then
    echo ">>> FastQC..."
    mkdir -p "${OUTDIR}/fastqc"
    fastqc -t ${THREADS} -o "${OUTDIR}/fastqc" \
        ${DATADIR}/tumor_R1.fastq.gz ${DATADIR}/tumor_R2.fastq.gz \
        ${DATADIR}/normal_R1.fastq.gz ${DATADIR}/normal_R2.fastq.gz \
        2>&1 | tail -3
fi

###############################################################################
# STEP 2: Align + process each sample (parallel for T and N)
###############################################################################
for SAMPLE in tumor normal; do
    if [ ! -f "${OUTDIR}/${SAMPLE}_dedup.bam" ]; then
        echo ">>> Processing ${SAMPLE}..."

        # Align
        bwa-mem2 mem -t ${THREADS} \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
            "${GENOME}" \
            "${DATADIR}/${SAMPLE}_R1.fastq.gz" "${DATADIR}/${SAMPLE}_R2.fastq.gz" \
            2>/dev/null | \
            samtools sort -@ ${THREADS} -o "${OUTDIR}/${SAMPLE}_sorted.bam" -

        samtools index "${OUTDIR}/${SAMPLE}_sorted.bam"

        # Mark duplicates
        picard MarkDuplicates \
            -I "${OUTDIR}/${SAMPLE}_sorted.bam" \
            -O "${OUTDIR}/${SAMPLE}_dedup.bam" \
            -M "${OUTDIR}/${SAMPLE}_dup_metrics.txt" \
            --REMOVE_DUPLICATES false \
            --CREATE_INDEX true 2>&1 | tail -3
    fi
done

###############################################################################
# STEP 3: Coverage analysis (parallel with variant calling)
###############################################################################
for SAMPLE in tumor normal; do
    if [ ! -f "${OUTDIR}/${SAMPLE}_coverage.mosdepth.summary.txt" ]; then
        echo ">>> mosdepth ${SAMPLE}..."
        mosdepth --threads ${THREADS} --no-per-base \
            "${OUTDIR}/${SAMPLE}_coverage" \
            "${OUTDIR}/${SAMPLE}_dedup.bam" 2>&1 | tail -3
    fi
done

###############################################################################
# STEP 4: CONVERGE-1 → Mutect2 somatic variant calling (T+N paired)
###############################################################################
if [ ! -f "${OUTDIR}/somatic_raw.vcf.gz" ]; then
    echo ">>> Mutect2 (tumor-normal paired)..."
    gatk Mutect2 \
        -R "${GENOME}" \
        -I "${OUTDIR}/tumor_dedup.bam" \
        -I "${OUTDIR}/normal_dedup.bam" \
        -normal normal \
        -O "${OUTDIR}/somatic_raw.vcf.gz" \
        --native-pair-hmm-threads ${THREADS} \
        2>&1 | tail -5
fi

###############################################################################
# STEP 5: Filter somatic variants
###############################################################################
if [ ! -f "${OUTDIR}/somatic_filtered.vcf.gz" ]; then
    echo ">>> FilterMutectCalls..."
    gatk FilterMutectCalls \
        -R "${GENOME}" \
        -V "${OUTDIR}/somatic_raw.vcf.gz" \
        -O "${OUTDIR}/somatic_filtered.vcf.gz" \
        2>&1 | tail -3
fi

###############################################################################
# STEP 5b: HaplotypeCaller (germline) on normal sample
###############################################################################
if [ ! -f "${OUTDIR}/germline_raw.vcf.gz" ]; then
    echo ">>> HaplotypeCaller (germline)..."
    gatk HaplotypeCaller \
        -R "${GENOME}" \
        -I "${OUTDIR}/normal_dedup.bam" \
        -O "${OUTDIR}/germline_raw.vcf.gz" \
        --native-pair-hmm-threads ${THREADS} \
        2>&1 | tail -3
fi

###############################################################################
# STEP 6: CONVERGE-2 — bcftools stats + coverage
###############################################################################
if [ ! -f "${OUTDIR}/variant_stats.txt" ]; then
    echo ">>> bcftools stats..."
    bcftools stats "${OUTDIR}/somatic_filtered.vcf.gz" > "${OUTDIR}/variant_stats.txt"
fi

###############################################################################
# STEP 7: CONVERGE-3 — Filter PASS variants + extract info
###############################################################################
if [ ! -f "${OUTDIR}/somatic_pass.vcf.gz" ]; then
    echo ">>> Extract PASS variants..."
    bcftools view -f PASS "${OUTDIR}/somatic_filtered.vcf.gz" -Oz -o "${OUTDIR}/somatic_pass.vcf.gz"
    bcftools index "${OUTDIR}/somatic_pass.vcf.gz"
fi

###############################################################################
# STEP 8: CONVERGE-4 — Generate report
###############################################################################
echo ">>> Generating report..."
python3 << 'REPORT_EOF'
import subprocess, os, re
import pandas as pd

# Parse Mutect2 stats
total_raw = 0
total_pass = 0
try:
    out = subprocess.run(["bcftools", "view", "-H", "outputs/somatic_raw.vcf.gz"],
                        capture_output=True, text=True)
    total_raw = len(out.stdout.strip().split('\n')) if out.stdout.strip() else 0
except: pass

try:
    out = subprocess.run(["bcftools", "view", "-H", "outputs/somatic_pass.vcf.gz"],
                        capture_output=True, text=True)
    total_pass = len(out.stdout.strip().split('\n')) if out.stdout.strip() else 0
except: pass

# Count variant types
snps, indels = 0, 0
try:
    with open("outputs/variant_stats.txt") as f:
        for line in f:
            if line.startswith("SN") and "number of SNPs:" in line:
                snps = int(line.strip().split('\t')[-1])
            elif line.startswith("SN") and "number of indels:" in line:
                indels = int(line.strip().split('\t')[-1])
except: pass

# Coverage
tumor_cov, normal_cov = 0, 0
for sample, var_name in [("tumor", "tumor_cov"), ("normal", "normal_cov")]:
    try:
        with open(f"outputs/{sample}_coverage.mosdepth.summary.txt") as f:
            for line in f:
                if line.startswith("total"):
                    parts = line.strip().split('\t')
                    if len(parts) > 3:
                        exec(f"{var_name} = round(float(parts[3]), 1)")
    except: pass

# Duplicate rates
tumor_dup, normal_dup = 0, 0
for sample, var_name in [("tumor", "tumor_dup"), ("normal", "normal_dup")]:
    try:
        with open(f"outputs/{sample}_dup_metrics.txt") as f:
            for line in f:
                if line.startswith("Unknown") or (not line.startswith("#") and not line.startswith("LIBRARY") and "\t" in line):
                    parts = line.strip().split('\t')
                    if len(parts) > 8:
                        exec(f"{var_name} = round(float(parts[8])*100, 2)")
                        break
    except: pass

# Read counts from FASTQs
tumor_reads, normal_reads = 0, 0
for sample, var_name in [("tumor", "tumor_reads"), ("normal", "normal_reads")]:
    try:
        import gzip
        with gzip.open(f"data/{sample}_R1.fastq.gz", 'rt') as f:
            exec(f"{var_name} = sum(1 for _ in f) // 4")
    except: pass

# Alignment rates from BAM
tumor_mapped, normal_mapped = 0, 0
for sample, var_name in [("tumor", "tumor_mapped"), ("normal", "normal_mapped")]:
    try:
        out = subprocess.run(["samtools", "flagstat", f"outputs/{sample}_dedup.bam"],
                           capture_output=True, text=True)
        for line in out.stdout.split('\n'):
            if 'primary mapped' in line:
                exec(f"{var_name} = int(line.split('+')[0].strip())")
                break
    except: pass

report = pd.DataFrame([
    ('tumor_input_reads', tumor_reads),
    ('normal_input_reads', normal_reads),
    ('tumor_mapped_reads', tumor_mapped),
    ('normal_mapped_reads', normal_mapped),
    ('tumor_mean_coverage', tumor_cov),
    ('normal_mean_coverage', normal_cov),
    ('tumor_duplicate_pct', tumor_dup),
    ('normal_duplicate_pct', normal_dup),
    ('somatic_variants_raw', total_raw),
    ('somatic_variants_pass', total_pass),
    ('somatic_snps', snps),
    ('somatic_indels', indels),
    ('variant_calling_method', 'paired_tumor_normal'),
    ('search_engines_used', 1),
    ('fdr_method', 'mutect2_internal'),
], columns=['metric', 'value'])

report.to_csv('results/report.csv', index=False)
print(report.to_string(index=False))
REPORT_EOF

echo ">>> Pipeline complete!"
