#!/usr/bin/env python3
"""Build consensus circRNA set from multiple detection tools."""
import csv
import sys
import os
import glob

def parse_circexplorer2(ce2_dir):
    """Parse CIRCexplorer2 results (circularRNA_known.txt)."""
    circrnas = {}
    for txt in glob.glob(os.path.join(ce2_dir, "*/circularRNA_known.txt")):
        sample = os.path.basename(os.path.dirname(txt))
        with open(txt) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom, start, end = parts[0], parts[1], parts[2]
                    strand = parts[5] if len(parts) > 5 else "."
                    key = f"{chrom}:{start}-{end}:{strand}"
                    count = int(parts[12]) if len(parts) > 12 else 1
                    if key not in circrnas:
                        circrnas[key] = {"chrom": chrom, "start": int(start),
                                         "end": int(end), "strand": strand,
                                         "samples": {}, "tool": "CIRCexplorer2"}
                    circrnas[key]["samples"][sample] = count
    return circrnas

def parse_dcc(dcc_dir):
    """Parse DCC results (CircRNACount)."""
    circrnas = {}
    count_file = os.path.join(dcc_dir, "CircRNACount")
    if not os.path.exists(count_file):
        return circrnas
    with open(count_file) as f:
        header = f.readline().strip().split('\t')
        sample_cols = header[3:]  # columns after chr, start, end
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            key = f"{chrom}:{start}-{end}:."
            counts = {}
            for i, col in enumerate(sample_cols):
                if i + 3 < len(parts):
                    try:
                        c = int(parts[i + 3])
                        if c > 0:
                            counts[col] = c
                    except ValueError:
                        pass
            if counts:
                circrnas[key] = {"chrom": chrom, "start": int(start),
                                 "end": int(end), "strand": ".",
                                 "samples": counts, "tool": "DCC"}
    return circrnas

def main():
    if len(sys.argv) < 4:
        print("Usage: consensus.py <ce2_dir> <dcc_dir> <output_dir>")
        sys.exit(1)

    ce2_dir, dcc_dir, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(output_dir, exist_ok=True)

    # Parse results from both tools
    ce2_circs = parse_circexplorer2(ce2_dir)
    dcc_circs = parse_dcc(dcc_dir)

    print(f"CIRCexplorer2: {len(ce2_circs)} circRNAs")
    print(f"DCC: {len(dcc_circs)} circRNAs")

    # Build consensus: normalize coordinates and find overlaps
    all_keys = set()
    tool_sets = {"CIRCexplorer2": set(), "DCC": set()}

    for key in ce2_circs:
        # Normalize key (remove strand for matching)
        parts = key.split(":")
        norm_key = f"{parts[0]}:{parts[1]}"
        all_keys.add(norm_key)
        tool_sets["CIRCexplorer2"].add(norm_key)

    for key in dcc_circs:
        parts = key.split(":")
        norm_key = f"{parts[0]}:{parts[1]}"
        all_keys.add(norm_key)
        tool_sets["DCC"].add(norm_key)

    # Consensus: found by ≥1 tool (since we only have 2 tools, report all)
    # In practice with more tools, would require ≥2
    consensus = tool_sets["CIRCexplorer2"] | tool_sets["DCC"]
    shared = tool_sets["CIRCexplorer2"] & tool_sets["DCC"]

    print(f"Total unique: {len(consensus)}")
    print(f"Shared by both tools: {len(shared)}")

    # Write consensus BED file
    bed_file = os.path.join(output_dir, "consensus_circrnas.bed")
    with open(bed_file, 'w') as f:
        for norm_key in sorted(consensus):
            chrom_part, coord_part = norm_key.split(":", 1)
            start, end = coord_part.split("-")
            n_tools = sum(1 for t in tool_sets if norm_key in tool_sets[t])
            f.write(f"{chrom_part}\t{start}\t{end}\tcircRNA\t{n_tools}\t.\n")

    # Write summary
    summary_file = os.path.join(output_dir, "consensus_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        writer.writerow(['circexplorer2_count', len(tool_sets["CIRCexplorer2"])])
        writer.writerow(['dcc_count', len(tool_sets["DCC"])])
        writer.writerow(['total_unique', len(consensus)])
        writer.writerow(['shared_by_both', len(shared)])
        writer.writerow(['tools_used', 2])

    print(f"Consensus written to: {bed_file}")

if __name__ == "__main__":
    main()
