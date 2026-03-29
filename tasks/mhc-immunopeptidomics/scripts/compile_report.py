#!/usr/bin/env python3
"""Compile MHC immunopeptidomics report."""
import csv, os, re
from collections import defaultdict

OUTDIR = "outputs"
SAMPLES = ["sample1", "sample2"]

def count_spectra(mzml_path):
    """Count spectra in mzML by counting <spectrum elements."""
    if not os.path.exists(mzml_path):
        return 0
    count = 0
    with open(mzml_path) as f:
        for line in f:
            if '<spectrum ' in line:
                count += 1
    return count

def parse_idxml_count(path):
    """Count peptide hits in idXML."""
    if not os.path.exists(path):
        return 0, set()
    count = 0
    peptides = set()
    with open(path) as f:
        for line in f:
            m = re.search(r'<PeptideHit.*?sequence="([^"]*)"', line)
            if m:
                count += 1
                peptides.add(m.group(1))
    return count, peptides

def parse_text_export(path):
    """Parse TextExporter TSV output."""
    if not os.path.exists(path):
        return []
    peptides = []
    with open(path) as f:
        in_peptide = False
        for line in f:
            if line.startswith('PEPTIDE'):
                in_peptide = True
                continue
            if in_peptide and line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    peptides.append(parts)
    return peptides

def count_features(featurexml_path):
    """Count features in featureXML."""
    if not os.path.exists(featurexml_path):
        return 0
    count = 0
    with open(featurexml_path) as f:
        for line in f:
            if '<feature id=' in line:
                count += 1
    return count

# Gather statistics
total_spectra = 0
total_psms_raw = 0
total_psms_filtered = 0
all_peptides_raw = set()
all_peptides_filtered = set()
total_features = 0
per_sample = {}

for sid in SAMPLES:
    # Spectra count
    spectra = count_spectra(f"data/{sid}.mzML")
    total_spectra += spectra

    # Raw search results
    raw_count, raw_peps = parse_idxml_count(f"{OUTDIR}/search/{sid}.idXML")
    total_psms_raw += raw_count
    all_peptides_raw.update(raw_peps)

    # Filtered results
    filt_count, filt_peps = parse_idxml_count(f"{OUTDIR}/filter/{sid}.idXML")
    total_psms_filtered += filt_count
    all_peptides_filtered.update(filt_peps)

    # Features
    feat_count = count_features(f"{OUTDIR}/quantify/{sid}.featureXML")
    total_features += feat_count

    per_sample[sid] = {
        'spectra': spectra,
        'psms_raw': raw_count,
        'psms_filtered': filt_count,
        'unique_peptides': len(filt_peps),
        'features': feat_count
    }

# Peptide length distribution (from filtered peptides)
length_dist = defaultdict(int)
for pep in all_peptides_filtered:
    length_dist[len(pep)] += 1

# Most common length
if length_dist:
    modal_length = max(length_dist, key=length_dist.get)
    pct_8_12 = sum(length_dist.get(l, 0) for l in range(8, 13))
    pct_8_12_frac = round(100 * pct_8_12 / len(all_peptides_filtered), 1) if all_peptides_filtered else 0
else:
    modal_length = 0
    pct_8_12_frac = 0

report = [
    ('total_samples', len(SAMPLES)),
    ('total_spectra', total_spectra),
    ('total_psms_searched', total_psms_raw),
    ('total_psms_filtered', total_psms_filtered),
    ('unique_peptides_raw', len(all_peptides_raw)),
    ('unique_peptides_filtered', len(all_peptides_filtered)),
    ('modal_peptide_length', modal_length),
    ('pct_peptides_8_12mer', pct_8_12_frac),
    ('total_quantified_features', total_features),
    ('sample1_spectra', per_sample['sample1']['spectra']),
    ('sample1_psms', per_sample['sample1']['psms_filtered']),
    ('sample2_spectra', per_sample['sample2']['spectra']),
    ('sample2_psms', per_sample['sample2']['psms_filtered']),
]

with open("results/report.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['metric', 'value'])
    writer.writerows(report)

print("=== Final Report ===")
for m, v in report:
    print(f"  {m} = {v}")
