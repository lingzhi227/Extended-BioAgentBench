#!/usr/bin/env python3
"""Merge multi-classifier results and generate report."""
import csv, os, json, math
from collections import defaultdict

OUTDIR = "outputs"
SAMPLES = ["sample1", "sample2"]

def parse_kreport(path):
    """Parse Kraken2/Centrifuge kreport format."""
    taxa = {}
    if not os.path.exists(path):
        return taxa
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                try:
                    pct = float(parts[0].strip())
                except ValueError:
                    continue
                rank = parts[3].strip()
                name = parts[5].strip()
                if rank in ('S', 'G') and pct > 0:
                    taxa[name] = {'pct': pct, 'rank': rank}
    return taxa

def parse_kaiju_table(path):
    """Parse kaiju2table output: file, percent, reads, taxon_id, taxon_name."""
    taxa = {}
    if not os.path.exists(path):
        return taxa
    with open(path) as f:
        header = True
        for line in f:
            if header:
                header = False
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                try:
                    pct = float(parts[1])
                except ValueError:
                    continue
                name = parts[4].strip()
                if name not in ('cannot be assigned', 'unclassified',
                                'cannot be assigned to a (non-viral) species',
                                'cannot be assigned to a (non-viral) genus') and pct > 0:
                    taxa[name] = {'pct': pct}
    return taxa

def parse_bracken(path):
    """Parse Bracken output."""
    taxa = {}
    if not os.path.exists(path):
        return taxa
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            name = row.get('name', '').strip()
            frac = float(row.get('fraction_total_reads', 0))
            if name and frac > 0:
                taxa[name] = {'pct': frac * 100}
    return taxa

def shannon(abundances):
    total = sum(abundances)
    if total == 0:
        return 0
    probs = [a/total for a in abundances if a > 0]
    return -sum(p * math.log(p) for p in probs)

def simpson(abundances):
    total = sum(abundances)
    if total == 0:
        return 0
    return 1 - sum((a/total)**2 for a in abundances if a > 0)

# ── Step 3: Merge classifier results ────────────────────────────────────────
print("[Step 3] Merging classifier results...")
for sid in SAMPLES:
    k2 = parse_kreport(f"{OUTDIR}/kraken2/{sid}.kreport")
    br = parse_bracken(f"{OUTDIR}/bracken/{sid}_S.bracken")
    cf = parse_kreport(f"{OUTDIR}/centrifuge/{sid}.kreport")
    kj = parse_kaiju_table(f"{OUTDIR}/kaiju/{sid}_species.tsv")

    all_taxa = set(list(k2.keys()) + list(br.keys()) + list(cf.keys()) + list(kj.keys()))

    merged = []
    for taxon in all_taxa:
        tools_detecting = 0
        abundances = []
        tool_list = []
        for name, d in [("classifier1", k2), ("bracken", br), ("classifier2", cf), ("classifier3", kj)]:
            if taxon in d:
                tools_detecting += 1
                abundances.append(d[taxon]['pct'])
                tool_list.append(name)
        mean_abund = sum(abundances) / len(abundances) if abundances else 0
        merged.append({
            'taxon': taxon,
            'tools_detected': tools_detecting,
            'mean_abundance': round(mean_abund, 4),
            'tools': ','.join(tool_list)
        })

    merged.sort(key=lambda x: -x['mean_abundance'])
    with open(f"{OUTDIR}/merged/{sid}_merged.tsv", 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['taxon', 'tools_detected', 'mean_abundance', 'tools'],
                               delimiter='\t')
        writer.writeheader()
        writer.writerows(merged)

    # Consensus: >=2 tools
    consensus = [m for m in merged if m['tools_detected'] >= 2]
    with open(f"{OUTDIR}/merged/{sid}_consensus.tsv", 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['taxon', 'tools_detected', 'mean_abundance', 'tools'],
                               delimiter='\t')
        writer.writeheader()
        writer.writerows(consensus)

    print(f"  {sid}: {len(merged)} total taxa, {len(consensus)} consensus (>=2 tools)")

# ── Step 4: Diversity ────────────────────────────────────────────────────────
print("[Step 4] Calculating diversity...")
div_results = []
for sid in SAMPLES:
    abunds = []
    with open(f"{OUTDIR}/merged/{sid}_merged.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            abunds.append(float(row['mean_abundance']))
    sh = shannon(abunds)
    si = simpson(abunds)
    rich = len(abunds)
    div_results.append({'sample': sid, 'shannon': round(sh, 4), 'simpson': round(si, 4), 'richness': rich})
    print(f"  {sid}: Shannon={sh:.4f}, Simpson={si:.4f}, Richness={rich}")

with open(f"{OUTDIR}/diversity/alpha_diversity.tsv", 'w') as f:
    writer = csv.DictWriter(f, fieldnames=['sample', 'shannon', 'simpson', 'richness'], delimiter='\t')
    writer.writeheader()
    writer.writerows(div_results)

# ── Step 5: Final report ────────────────────────────────────────────────────
print("[Step 5] Generating report...")

total_input = 0
total_qc = 0
for sid in SAMPLES:
    with open(f"{OUTDIR}/fastp/{sid}_fastp.json") as f:
        j = json.load(f)
        total_input += j['summary']['before_filtering']['total_reads']
        total_qc += j['summary']['after_filtering']['total_reads']

# Count classified reads per tool from kreports
def count_kreport_classified(path):
    """Count classified reads from kreport. Returns (classified, unclassified)."""
    if not os.path.exists(path):
        return 0, 0
    unclassified = 0
    total = 0
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                try:
                    count = int(parts[1].strip())
                except ValueError:
                    continue
                rank = parts[3].strip()
                if rank == 'U':
                    unclassified = count
                    total += count
                elif rank in ('R', '-') and parts[5].strip() == 'root':
                    total += count
    classified = total - unclassified
    return classified, unclassified

k2_class = sum(count_kreport_classified(f"{OUTDIR}/kraken2/{s}.kreport")[0] for s in SAMPLES)
cf_class = sum(count_kreport_classified(f"{OUTDIR}/centrifuge/{s}.kreport")[0] for s in SAMPLES)

# Kaiju classified reads
kj_class = 0
for sid in SAMPLES:
    path = f"{OUTDIR}/kaiju/{sid}.out"
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                if line.startswith('C'):
                    kj_class += 1

# Species counts per tool
k2_species = set()
cf_species = set()
kj_species = set()
for sid in SAMPLES:
    for taxon in parse_kreport(f"{OUTDIR}/kraken2/{sid}.kreport"):
        k2_species.add(taxon)
    for taxon in parse_kreport(f"{OUTDIR}/centrifuge/{sid}.kreport"):
        cf_species.add(taxon)
    for taxon in parse_kaiju_table(f"{OUTDIR}/kaiju/{sid}_species.tsv"):
        kj_species.add(taxon)

# Total unique and consensus
all_species = set()
consensus_species = set()
for sid in SAMPLES:
    with open(f"{OUTDIR}/merged/{sid}_merged.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            all_species.add(row['taxon'])
    with open(f"{OUTDIR}/merged/{sid}_consensus.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            consensus_species.add(row['taxon'])

# Top taxon
top_taxon = ""
top_abund = 0
for sid in SAMPLES:
    with open(f"{OUTDIR}/merged/{sid}_merged.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            a = float(row['mean_abundance'])
            if a > top_abund:
                top_abund = a
                top_taxon = row['taxon']

mean_shannon = round(sum(d['shannon'] for d in div_results) / len(div_results), 4) if div_results else 0
mean_simpson = round(sum(d['simpson'] for d in div_results) / len(div_results), 4) if div_results else 0
mean_richness = round(sum(d['richness'] for d in div_results) / len(div_results), 1) if div_results else 0

report = [
    ('total_samples', len(SAMPLES)),
    ('total_input_reads', total_input),
    ('total_qc_reads', total_qc),
    ('classifiers_used', 3),
    ('classified_reads_classifier1', k2_class),
    ('classified_reads_classifier2', cf_class),
    ('classified_reads_classifier3', kj_class),
    ('unique_species_classifier1', len(k2_species)),
    ('unique_species_classifier2', len(cf_species)),
    ('unique_species_classifier3', len(kj_species)),
    ('total_unique_taxa', len(all_species)),
    ('consensus_taxa', len(consensus_species)),
    ('top_taxon', top_taxon if top_taxon else 'none'),
    ('top_taxon_abundance', round(top_abund, 2)),
    ('mean_shannon', mean_shannon),
    ('mean_simpson', mean_simpson),
    ('mean_richness', mean_richness),
]

with open("results/report.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['metric', 'value'])
    writer.writerows(report)

print("=== Final Report ===")
for m, v in report:
    print(f"  {m} = {v}")
