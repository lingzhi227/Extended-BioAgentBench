#!/usr/bin/env python3
"""Spectral library matching using matchms."""
import csv
import sys
import os
import json

def parse_msp(msp_file):
    """Parse MSP spectral library file."""
    spectra = []
    current = {}
    peaks_mz = []
    peaks_int = []
    in_peaks = False

    with open(msp_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                if current and peaks_mz:
                    current['mz'] = peaks_mz
                    current['intensities'] = peaks_int
                    spectra.append(current)
                current = {}
                peaks_mz = []
                peaks_int = []
                in_peaks = False
                continue

            if line.lower().startswith('num peaks'):
                in_peaks = True
                continue
            elif ':' in line and not in_peaks:
                key, val = line.split(':', 1)
                current[key.strip()] = val.strip()
            elif in_peaks:
                # Parse peak pairs: "mz intensity; mz intensity; ..."
                pairs = line.split(';')
                for pair in pairs:
                    pair = pair.strip()
                    if pair:
                        parts = pair.split()
                        if len(parts) >= 2:
                            try:
                                peaks_mz.append(float(parts[0]))
                                peaks_int.append(float(parts[1]))
                            except ValueError:
                                pass

    # Last spectrum
    if current and peaks_mz:
        current['mz'] = peaks_mz
        current['intensities'] = peaks_int
        spectra.append(current)

    return spectra

def cosine_score(mz1, int1, mz2, int2, tolerance=0.5):
    """Compute cosine similarity between two spectra."""
    import math

    if not mz1 or not mz2:
        return 0.0, 0

    # Normalize intensities
    max1 = max(int1) if int1 else 1
    max2 = max(int2) if int2 else 1
    norm1 = [i / max1 for i in int1]
    norm2 = [i / max2 for i in int2]

    # Find matching peaks
    matched_prod = 0.0
    n_matches = 0
    used2 = set()

    for i, m1 in enumerate(mz1):
        best_j = -1
        best_diff = tolerance + 1
        for j, m2 in enumerate(mz2):
            if j in used2:
                continue
            diff = abs(m1 - m2)
            if diff <= tolerance and diff < best_diff:
                best_diff = diff
                best_j = j
        if best_j >= 0:
            matched_prod += norm1[i] * norm2[best_j]
            n_matches += 1
            used2.add(best_j)

    # Cosine similarity
    sum_sq1 = sum(x**2 for x in norm1)
    sum_sq2 = sum(x**2 for x in norm2)

    if sum_sq1 == 0 or sum_sq2 == 0:
        return 0.0, 0

    score = matched_prod / (math.sqrt(sum_sq1) * math.sqrt(sum_sq2))
    return round(score, 4), n_matches

def load_pseudospectra(ps_file, peaklist_file):
    """Load pseudospectra from CAMERA output."""
    # Load pseudospectra info
    ps_info = []
    with open(ps_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ps_info.append(row)

    # Group by pcgroup
    groups = {}
    for row in ps_info:
        grp = row['pcgroup']
        if grp not in groups:
            groups[grp] = {'mz': [], 'intensities': []}
        groups[grp]['mz'].append(float(row['mz']))

    # Load intensities from peaklist
    with open(peaklist_file) as f:
        reader = csv.DictReader(f)
        peaklist = list(reader)

    # Add mean intensities to groups
    for i, row in enumerate(ps_info):
        grp = row['pcgroup']
        # Use mean intensity across samples as proxy
        pk = peaklist[i] if i < len(peaklist) else {}
        # Find numeric columns (sample intensities)
        intensity = 0
        n_vals = 0
        for k, v in pk.items():
            try:
                val = float(v)
                if val > 0 and k not in ['mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax',
                                           'npeaks', 'pcgroup']:
                    intensity += val
                    n_vals += 1
            except (ValueError, TypeError):
                pass
        if n_vals > 0:
            groups[grp]['intensities'].append(intensity / n_vals)
        else:
            groups[grp]['intensities'].append(0)

    return groups

def main():
    if len(sys.argv) < 5:
        print("Usage: matchms_match.py <library.msp> <pseudospectra_info.csv> <camera_peaklist.csv> <output_dir>")
        sys.exit(1)

    lib_file = sys.argv[1]
    ps_file = sys.argv[2]
    peaklist_file = sys.argv[3]
    output_dir = sys.argv[4]
    os.makedirs(output_dir, exist_ok=True)

    # Load library
    library = parse_msp(lib_file)
    print(f"Loaded {len(library)} library spectra")
    for s in library:
        print(f"  - {s.get('Name', 'Unknown')}: {len(s.get('mz', []))} peaks")

    # Load pseudospectra
    groups = load_pseudospectra(ps_file, peaklist_file)
    print(f"Loaded {len(groups)} pseudospectra groups")

    # Match each pseudospectrum against library
    matches = []
    for grp_id, grp_data in sorted(groups.items(), key=lambda x: int(x[0])):
        best_score = 0
        best_match = ""
        best_n_matches = 0

        for lib_spec in library:
            score, n_matches = cosine_score(
                grp_data['mz'], grp_data['intensities'],
                lib_spec['mz'], lib_spec['intensities'],
                tolerance=0.5
            )
            if score > best_score:
                best_score = score
                best_match = lib_spec.get('Name', 'Unknown')
                best_n_matches = n_matches

        if best_score >= 0.3:  # minimum match threshold
            matches.append({
                'pseudospectrum_group': grp_id,
                'matched_compound': best_match,
                'cosine_score': best_score,
                'matched_peaks': best_n_matches,
                'query_peaks': len(grp_data['mz'])
            })

    # Write matches
    out_file = os.path.join(output_dir, "library_matches.csv")
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['pseudospectrum_group', 'matched_compound',
                                                'cosine_score', 'matched_peaks', 'query_peaks'])
        writer.writeheader()
        writer.writerows(matches)

    # Summary
    n_matched = len(matches)
    compounds = list(set(m['matched_compound'] for m in matches))
    mean_score = sum(m['cosine_score'] for m in matches) / len(matches) if matches else 0

    summary_file = os.path.join(output_dir, "matching_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        writer.writerow(['library_spectra', len(library)])
        writer.writerow(['pseudospectra_matched', n_matched])
        writer.writerow(['unique_compounds_identified', len(compounds)])
        writer.writerow(['mean_cosine_score', f"{mean_score:.3f}"])
        writer.writerow(['identified_compounds', ';'.join(sorted(compounds))])

    print(f"\n=== Matching Results ===")
    print(f"Matched: {n_matched}/{len(groups)} pseudospectra")
    print(f"Unique compounds: {len(compounds)}")
    if compounds:
        for c in sorted(compounds):
            print(f"  - {c}")

if __name__ == "__main__":
    main()
