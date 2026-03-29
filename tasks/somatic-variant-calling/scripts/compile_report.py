#!/usr/bin/env python3
"""Compile somatic variant calling report from pipeline outputs."""
import csv, os, json, subprocess

OUTDIR = "outputs"

def count_vcf_variants(path):
    """Count variants in a VCF file (excluding header lines)."""
    if not os.path.exists(path):
        return 0
    count = 0
    import gzip
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    count += 1
    except Exception:
        return 0
    return count

def count_pass_variants(path):
    """Count PASS variants in a VCF."""
    if not os.path.exists(path):
        return 0
    count = 0
    import gzip
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.split('\t')
                    if len(fields) >= 7 and (fields[6] == 'PASS' or fields[6] == '.'):
                        count += 1
    except Exception:
        return 0
    return count

def parse_mosdepth_summary(path):
    """Parse mosdepth summary file."""
    if not os.path.exists(path):
        return {'mean': 0, 'min': 0, 'max': 0}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('chrom') == 'total' or row.get('chrom', '').startswith('total'):
                return {
                    'mean': float(row.get('mean', 0)),
                    'min': float(row.get('min', 0)),
                    'max': float(row.get('max', 0))
                }
    return {'mean': 0, 'min': 0, 'max': 0}

def parse_fastp_json(path):
    """Parse fastp JSON report."""
    if not os.path.exists(path):
        return {'before': 0, 'after': 0}
    with open(path) as f:
        j = json.load(f)
    return {
        'before': j['summary']['before_filtering']['total_reads'],
        'after': j['summary']['after_filtering']['total_reads']
    }

def parse_markdup_metrics(path):
    """Parse GATK MarkDuplicates metrics."""
    if not os.path.exists(path):
        return {'dup_pct': 0}
    with open(path) as f:
        in_metrics = False
        headers = []
        for line in f:
            if line.startswith('## METRICS CLASS'):
                in_metrics = True
                continue
            if in_metrics and not headers:
                headers = line.strip().split('\t')
                continue
            if in_metrics and headers:
                vals = line.strip().split('\t')
                if len(vals) >= len(headers):
                    d = dict(zip(headers, vals))
                    return {'dup_pct': float(d.get('PERCENT_DUPLICATION', 0))}
                break
    return {'dup_pct': 0}

# Gather results
tumor_fastp = parse_fastp_json(f"{OUTDIR}/fastp/tumor_fastp.json")
normal_fastp = parse_fastp_json(f"{OUTDIR}/fastp/normal_fastp.json")

tumor_dup = parse_markdup_metrics(f"{OUTDIR}/markdup/tumor.metrics.txt")
normal_dup = parse_markdup_metrics(f"{OUTDIR}/markdup/normal.metrics.txt")

tumor_cov = parse_mosdepth_summary(f"{OUTDIR}/coverage/tumor.mosdepth.summary.txt")
normal_cov = parse_mosdepth_summary(f"{OUTDIR}/coverage/normal.mosdepth.summary.txt")

# Variant counts per caller
caller1_raw = count_vcf_variants(f"{OUTDIR}/mutect2/somatic_raw.vcf.gz")
caller1_pass = count_pass_variants(f"{OUTDIR}/mutect2/somatic_filtered.vcf.gz")
caller2_raw = count_vcf_variants(f"{OUTDIR}/freebayes/joint_raw.vcf")
caller3_raw = count_vcf_variants(f"{OUTDIR}/bcftools_call/pileup_raw.vcf.gz")

# Annotated counts
ann1 = count_vcf_variants(f"{OUTDIR}/annotate/mutect2_annotated.vcf.gz")
ann2 = count_vcf_variants(f"{OUTDIR}/annotate/freebayes_annotated.vcf.gz")
ann3 = count_vcf_variants(f"{OUTDIR}/annotate/bcftools_annotated.vcf.gz")

# Contamination
contamination = 0.0
contam_file = f"{OUTDIR}/mutect2/contamination.table"
if os.path.exists(contam_file):
    with open(contam_file) as f:
        for line in f:
            if not line.startswith('sample'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    contamination = float(parts[1])

report = [
    ('tumor_input_reads', tumor_fastp['before']),
    ('tumor_qc_reads', tumor_fastp['after']),
    ('normal_input_reads', normal_fastp['before']),
    ('normal_qc_reads', normal_fastp['after']),
    ('tumor_duplication_rate', round(tumor_dup['dup_pct'], 4)),
    ('normal_duplication_rate', round(normal_dup['dup_pct'], 4)),
    ('tumor_mean_coverage', round(tumor_cov['mean'], 2)),
    ('normal_mean_coverage', round(normal_cov['mean'], 2)),
    ('caller1_raw_variants', caller1_raw),
    ('caller1_pass_variants', caller1_pass),
    ('caller2_raw_variants', caller2_raw),
    ('caller3_raw_variants', caller3_raw),
    ('callers_used', 3),
    ('annotated_variants_caller1', ann1),
    ('annotated_variants_caller2', ann2),
    ('annotated_variants_caller3', ann3),
    ('estimated_contamination', round(contamination, 4)),
]

with open("results/report.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['metric', 'value'])
    writer.writerows(report)

print("=== Final Report ===")
for m, v in report:
    print(f"  {m} = {v}")
