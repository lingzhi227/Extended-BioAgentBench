#!/usr/bin/env python3
"""Parse pairtools dedup stats for contact statistics."""
import os, json

OUTDIR = "outputs"
SAMPLES = ["KD_rep1", "KD_rep2", "WT_rep1", "WT_rep2"]

stats = {}
for sid in SAMPLES:
    stats_file = f"{OUTDIR}/pairs/{sid}.dedup.stats"
    s = {'total': 0, 'cis': 0, 'trans': 0, 'cis_1kb': 0, 'dups': 0}
    if os.path.exists(stats_file):
        with open(stats_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t') if '\t' in line else line.split()
                if len(parts) >= 2:
                    key, val = parts[0], parts[1]
                    try:
                        val = int(val)
                    except ValueError:
                        continue
                    if 'total' in key.lower():
                        s['total'] = val
                    elif key.lower() in ('cis', 'cis_1kb+', 'cis_10kb+'):
                        s['cis'] = max(s['cis'], val)
                    elif key.lower() == 'trans':
                        s['trans'] = val
                    elif 'dup' in key.lower():
                        s['dups'] = val
    stats[sid] = s

with open(f"{OUTDIR}/stats/contact_stats.json", 'w') as f:
    json.dump(stats, f, indent=2)

for sid, s in stats.items():
    print(f"  {sid}: total={s['total']}, cis={s['cis']}, trans={s['trans']}")
