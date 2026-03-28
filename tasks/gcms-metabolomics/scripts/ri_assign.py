#!/usr/bin/env python3
"""Compute Kovats retention indices from alkane reference standards."""
import csv
import sys
import math
import os

def load_alkane_references(ref_file):
    """Load alkane RT-RI reference pairs."""
    refs = []
    with open(ref_file) as f:
        # Space-separated: rt RI (with header)
        lines = f.readlines()
        for line in lines[1:]:  # skip header
            parts = line.strip().split()
            if len(parts) >= 3:
                # Format: index rt RI
                rt = float(parts[1])
                ri = float(parts[2])
                refs.append((rt, ri))
    refs.sort(key=lambda x: x[0])
    return refs

def compute_kovats_ri(rt, refs):
    """Compute Kovats retention index for a given RT using linear interpolation."""
    if rt <= refs[0][0]:
        return refs[0][1]
    if rt >= refs[-1][0]:
        return refs[-1][1]
    for i in range(len(refs) - 1):
        rt_n, ri_n = refs[i]
        rt_n1, ri_n1 = refs[i + 1]
        if rt_n <= rt <= rt_n1:
            if rt_n1 == rt_n:
                return ri_n
            ri = ri_n + (ri_n1 - ri_n) * (rt - rt_n) / (rt_n1 - rt_n)
            return round(ri, 1)
    return None

def main():
    if len(sys.argv) < 4:
        print("Usage: ri_assign.py <feature_defs.csv> <alkane_ref.csv> <output_dir>")
        sys.exit(1)

    feat_file = sys.argv[1]
    ref_file = sys.argv[2]
    output_dir = sys.argv[3]
    os.makedirs(output_dir, exist_ok=True)

    # Load alkane references
    refs = load_alkane_references(ref_file)
    print(f"Loaded {len(refs)} alkane reference points")
    print(f"RT range: {refs[0][0]:.2f} - {refs[-1][0]:.2f} min")
    print(f"RI range: {refs[0][1]:.0f} - {refs[-1][1]:.0f}")

    # Load feature definitions (from XCMS)
    features = []
    with open(feat_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            features.append(row)

    print(f"Processing {len(features)} features")

    # Compute RI for each feature using median RT
    ri_results = []
    ri_values = []
    for feat in features:
        feat_id = feat.get('', feat.get('X', 'unknown'))
        # RT is in seconds in XCMS, convert to minutes
        rt_sec = float(feat.get('rtmed', feat.get('rt', 0)))
        rt_min = rt_sec / 60.0

        ri = compute_kovats_ri(rt_min, refs)
        ri_results.append({
            'feature_id': feat_id,
            'rt_seconds': round(rt_sec, 2),
            'rt_minutes': round(rt_min, 2),
            'retention_index': ri
        })
        if ri is not None:
            ri_values.append(ri)

    # Write RI assignments
    out_file = os.path.join(output_dir, "ri_assignments.csv")
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['feature_id', 'rt_seconds', 'rt_minutes', 'retention_index'])
        writer.writeheader()
        writer.writerows(ri_results)

    # Summary
    n_assigned = len([r for r in ri_results if r['retention_index'] is not None])
    ri_min = min(ri_values) if ri_values else 0
    ri_max = max(ri_values) if ri_values else 0
    ri_mean = sum(ri_values) / len(ri_values) if ri_values else 0

    summary_file = os.path.join(output_dir, "ri_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        writer.writerow(['features_with_ri', n_assigned])
        writer.writerow(['ri_range', f"{ri_min:.0f}-{ri_max:.0f}"])
        writer.writerow(['mean_ri', f"{ri_mean:.1f}"])

    print(f"RI assigned to {n_assigned}/{len(ri_results)} features")
    print(f"RI range: {ri_min:.0f} - {ri_max:.0f}")
    print(f"Output: {out_file}")

if __name__ == "__main__":
    main()
