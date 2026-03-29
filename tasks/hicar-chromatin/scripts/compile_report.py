#!/usr/bin/env python3
"""Compile HiCAR chromatin interaction report."""
import csv, os, json

OUTDIR = "outputs"
SAMPLES = ["KD_rep1", "KD_rep2", "WT_rep1", "WT_rep2"]

# Read contact stats
contact_stats = {}
cs_file = f"{OUTDIR}/stats/contact_stats.json"
if os.path.exists(cs_file):
    with open(cs_file) as f:
        contact_stats = json.load(f)

# Count total reads from cutadapt logs
total_input = 0
total_trimmed = 0
for sid in SAMPLES:
    log = f"{OUTDIR}/trimmed/{sid}.log"
    if os.path.exists(log):
        with open(log) as f:
            for line in f:
                if 'Total read pairs processed' in line:
                    val = line.split(':')[-1].strip().replace(',', '')
                    total_input += int(val)
                elif 'Pairs written (passing filters)' in line:
                    val = line.split(':')[-1].strip().split()[0].replace(',', '')
                    total_trimmed += int(val)

# Pair statistics
total_pairs = sum(s.get('total', 0) for s in contact_stats.values())
total_cis = sum(s.get('cis', 0) for s in contact_stats.values())
total_trans = sum(s.get('trans', 0) for s in contact_stats.values())
cis_ratio = round(total_cis / max(total_pairs, 1) * 100, 2)

# Peak counts
total_peaks = 0
per_sample_peaks = {}
for sid in SAMPLES:
    pf = f"{OUTDIR}/stats/{sid}_peak_count.txt"
    if os.path.exists(pf):
        with open(pf) as f:
            count = int(f.read().strip() or 0)
    else:
        count = 0
    per_sample_peaks[sid] = count
    total_peaks += count

# Cooler files
matrices_created = sum(1 for sid in SAMPLES
                       if os.path.exists(f"{OUTDIR}/cooler/{sid}.cool")
                       and os.path.getsize(f"{OUTDIR}/cooler/{sid}.cool") > 0)

report = [
    ('total_samples', len(SAMPLES)),
    ('conditions', 2),
    ('total_read_pairs', total_input),
    ('trimmed_read_pairs', total_trimmed),
    ('total_valid_pairs', total_pairs),
    ('cis_contacts', total_cis),
    ('trans_contacts', total_trans),
    ('cis_contact_pct', cis_ratio),
    ('total_accessibility_peaks', total_peaks),
    ('kd_rep1_peaks', per_sample_peaks.get('KD_rep1', 0)),
    ('kd_rep2_peaks', per_sample_peaks.get('KD_rep2', 0)),
    ('wt_rep1_peaks', per_sample_peaks.get('WT_rep1', 0)),
    ('wt_rep2_peaks', per_sample_peaks.get('WT_rep2', 0)),
    ('contact_matrices_created', matrices_created),
    ('matrix_resolution', 10000),
]

with open("results/report.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['metric', 'value'])
    writer.writerows(report)

print("=== Final Report ===")
for m, v in report:
    print(f"  {m} = {v}")
