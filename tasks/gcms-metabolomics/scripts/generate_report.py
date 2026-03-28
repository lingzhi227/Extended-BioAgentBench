#!/usr/bin/env python3
"""Merge all metabolomics results into final report CSV."""
import csv
import sys
import os
import math

def read_csv_dict(filepath):
    """Read a metric,value CSV into a dict."""
    result = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            result[row['metric']] = row['value']
    return result

def compute_sample_correlation(peak_matrix_file):
    """Compute mean pairwise correlation between samples."""
    # Read peak matrix
    with open(peak_matrix_file) as f:
        reader = csv.reader(f)
        header = next(reader)
        # First column is feature ID
        sample_names = header[1:]
        n_samples = len(sample_names)

        # Collect intensities per sample
        sample_data = {s: [] for s in sample_names}
        for row in reader:
            for i, s in enumerate(sample_names):
                try:
                    val = float(row[i + 1]) if row[i + 1] and row[i + 1] != 'NA' else 0
                except (ValueError, IndexError):
                    val = 0
                sample_data[s].append(val)

    # Compute pairwise Pearson correlations
    correlations = []
    samples = list(sample_data.keys())
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            x = sample_data[samples[i]]
            y = sample_data[samples[j]]
            n = len(x)
            if n < 2:
                continue
            mean_x = sum(x) / n
            mean_y = sum(y) / n
            cov = sum((x[k] - mean_x) * (y[k] - mean_y) for k in range(n))
            var_x = sum((x[k] - mean_x) ** 2 for k in range(n))
            var_y = sum((y[k] - mean_y) ** 2 for k in range(n))
            if var_x > 0 and var_y > 0:
                r = cov / (math.sqrt(var_x) * math.sqrt(var_y))
                correlations.append(r)

    mean_corr = sum(correlations) / len(correlations) if correlations else 0
    return round(mean_corr, 3)

def compute_cv(peak_matrix_file):
    """Compute mean CV across features."""
    with open(peak_matrix_file) as f:
        reader = csv.reader(f)
        header = next(reader)
        n_samples = len(header) - 1

        cvs = []
        for row in reader:
            vals = []
            for v in row[1:]:
                try:
                    val = float(v) if v and v != 'NA' else None
                except ValueError:
                    val = None
                if val is not None and val > 0:
                    vals.append(val)
            if len(vals) >= 2:
                mean_v = sum(vals) / len(vals)
                if mean_v > 0:
                    std_v = math.sqrt(sum((v - mean_v)**2 for v in vals) / (len(vals) - 1))
                    cvs.append(std_v / mean_v * 100)

    mean_cv = sum(cvs) / len(cvs) if cvs else 0
    return round(mean_cv, 1)

def main():
    if len(sys.argv) < 3:
        print("Usage: generate_report.py <outputs_dir> <results_dir>")
        sys.exit(1)

    outputs_dir = sys.argv[1]
    results_dir = sys.argv[2]
    os.makedirs(results_dir, exist_ok=True)

    # Read all summaries
    xcms_sum = read_csv_dict(os.path.join(outputs_dir, "xcms", "xcms_summary.csv"))
    camera_sum = read_csv_dict(os.path.join(outputs_dir, "camera", "camera_summary.csv"))
    ri_sum = read_csv_dict(os.path.join(outputs_dir, "ri", "ri_summary.csv"))
    match_sum = read_csv_dict(os.path.join(outputs_dir, "matching", "matching_summary.csv"))

    # Compute additional stats
    peak_matrix = os.path.join(outputs_dir, "xcms", "peak_matrix.csv")
    mean_corr = compute_sample_correlation(peak_matrix)
    mean_cv = compute_cv(peak_matrix)

    # Build final report
    report = [
        ("total_peaks_detected", xcms_sum.get("total_peaks_detected", "")),
        ("features_before_alignment", xcms_sum.get("features_before_alignment", "")),
        ("features_after_alignment", xcms_sum.get("features_after_alignment", "")),
        ("samples_processed", xcms_sum.get("samples_processed", "")),
        ("fill_rate_percent", xcms_sum.get("fill_rate_percent", "")),
        ("pseudospectra_count", camera_sum.get("pseudospectra_count", "")),
        ("isotope_annotations", camera_sum.get("isotope_annotations", "")),
        ("adduct_annotations", camera_sum.get("adduct_annotations", "")),
        ("retention_index_range", ri_sum.get("ri_range", "")),
        ("mean_retention_index", ri_sum.get("mean_ri", "")),
        ("library_spectra_count", match_sum.get("library_spectra", "")),
        ("compounds_matched", match_sum.get("pseudospectra_matched", "")),
        ("unique_compounds_identified", match_sum.get("unique_compounds_identified", "")),
        ("mean_match_score", match_sum.get("mean_cosine_score", "")),
        ("identified_compounds", match_sum.get("identified_compounds", "")),
        ("mean_sample_correlation", str(mean_corr)),
        ("mean_feature_cv_percent", str(mean_cv)),
    ]

    # Write report
    out_file = os.path.join(results_dir, "report.csv")
    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        for metric, value in report:
            writer.writerow([metric, value])

    print("=== Final Report ===")
    for metric, value in report:
        print(f"  {metric}: {value}")
    print(f"\nReport written to: {out_file}")

if __name__ == "__main__":
    main()
