#!/usr/bin/env python3
"""Quantify circRNAs across samples and generate report."""
import csv
import sys
import os
import glob

def main():
    if len(sys.argv) < 5:
        print("Usage: quantify_and_report.py <consensus.bed> <ce2_dir> <genome.fa> <results_dir>")
        sys.exit(1)

    bed_file, ce2_dir, genome_fa, results_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    os.makedirs(results_dir, exist_ok=True)

    # Read consensus circRNAs
    circrnas = []
    with open(bed_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                circrnas.append({
                    'chrom': parts[0], 'start': int(parts[1]),
                    'end': int(parts[2]), 'n_tools': int(parts[4])
                })

    # Collect per-sample counts from CIRCexplorer2
    sample_counts = {}
    for txt in sorted(glob.glob(os.path.join(ce2_dir, "*/circularRNA_known.txt"))):
        sample = os.path.basename(os.path.dirname(txt))
        sample_counts[sample] = {}
        with open(txt) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 13:
                    key = f"{parts[0]}:{parts[1]}-{parts[2]}"
                    count = int(parts[12])
                    sample_counts[sample][key] = count

    # Compute statistics
    n_consensus = len(circrnas)
    n_shared = sum(1 for c in circrnas if c['n_tools'] >= 2)

    # Length distribution
    lengths = [c['end'] - c['start'] for c in circrnas]
    mean_length = sum(lengths) / len(lengths) if lengths else 0
    min_length = min(lengths) if lengths else 0
    max_length = max(lengths) if lengths else 0

    # Count total BSJ reads
    total_bsj = 0
    for sample, counts in sample_counts.items():
        total_bsj += sum(counts.values())

    # Samples processed
    n_samples = len(sample_counts)

    # Extract sequences using bedtools (compute length stats only)
    # Chromosome distribution
    chroms = {}
    for c in circrnas:
        chroms[c['chrom']] = chroms.get(c['chrom'], 0) + 1

    # Write report
    report_file = os.path.join(results_dir, "report.csv")
    with open(report_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        writer.writerow(['samples_analyzed', n_samples])
        writer.writerow(['total_circular_rnas', n_consensus])
        writer.writerow(['shared_by_multiple_tools', n_shared])
        writer.writerow(['total_back_splice_reads', total_bsj])
        writer.writerow(['mean_circrna_length', int(mean_length)])
        writer.writerow(['min_circrna_length', min_length])
        writer.writerow(['max_circrna_length', max_length])
        writer.writerow(['detection_tools_used', 2])
        writer.writerow(['chromosomes_with_circrnas', len(chroms)])
        # Per-condition summary
        n2_count = sum(1 for s in sample_counts if 'N2' in s)
        fust1_count = sum(1 for s in sample_counts if 'fust1' in s)
        writer.writerow(['n2_samples', n2_count])
        writer.writerow(['fust1_samples', fust1_count])

    print(f"Report written to: {report_file}")
    with open(report_file) as f:
        print(f.read())

if __name__ == "__main__":
    main()
