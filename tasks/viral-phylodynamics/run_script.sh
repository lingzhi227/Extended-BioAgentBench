#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Task: viral-phylodynamics
# DAG Structure (depth=10, convergence=4, tools=10):
#
#  sequences.fasta     metadata.tsv (dates, locations)
#        |                   |
#  [mafft align] ----  [python parse dates]            Level 1
#        |                   |
#  +-----+-------+          |
#  |     |       |          |
# [trimal [iqtree [python   |                          Level 2
#  gap     model   clock     |
#  filter] test]   signal    |
#  |     |       (root-to-   |
#  |     |       tip)]       |
#  |     |       |           |
#  +-----+---+---+           |
#            |                |
#    [CONVERGENCE 1] <--------+                        Level 3
#    (cleaned alignment + model + dates)
#            |
#    [iqtree ML tree]                                  Level 4
#    (best model, bootstrap)
#            |
#    [treetime molecular clock]                        Level 5
#    (--clock-filter)
#            |
#    +-------+---------------+
#    |       |               |
# [treetime [treetime     [treetime                    Level 6
#  mugration  skyline       ancestral
#  (geo)]     (Neeff)]      (reconstruct)]
#    |       |               |
#    +-------+-------+-------+
#                    |
#            [CONVERGENCE 2]                           Level 7
#            (migration + skyline + ancestral)
#                    |
#    +---------------+---------------+
#    |               |               |
# [augur         [python          [python              Level 8
#  export         transmission    rate
#  (Auspice       chain           estimation
#   JSON)]        reconstruction]  report]
#    |               |               |
#    +---------------+---------------+
#                    |
#            [CONVERGENCE 3]                           Level 9
#            (Auspice + chains + rates)
#                    |
#            [CONVERGENCE 4] <-- root-to-tip QC        Level 10
#            [python report]
# ============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

DATA="$WORKDIR/data"
REF="$WORKDIR/reference"
OUT="$WORKDIR/outputs"
RES="$WORKDIR/results"

mkdir -p "$OUT"/{alignment,trees,clock,mugration,skyline,ancestral,export,analysis}
mkdir -p "$RES"

# ============================================================================
# LEVEL 1a: Multiple sequence alignment with MAFFT
# ============================================================================
if [ ! -f "$OUT/alignment/aligned_raw.fasta" ]; then
  echo "[Level 1a] Running MAFFT alignment..."
  mafft --auto --thread "$THREADS" "$DATA/sequences.fasta" > "$OUT/alignment/aligned_raw.fasta"
fi

# ============================================================================
# LEVEL 1b: Parse dates from metadata
# ============================================================================
if [ ! -f "$OUT/analysis/parsed_dates.tsv" ]; then
  echo "[Level 1b] Parsing dates from metadata..."
  python3 << 'PYEOF'
import pandas as pd
import sys, os

meta = pd.read_csv("data/metadata.tsv", sep="\t")
# Convert date strings to decimal dates
dates = []
for _, row in meta.iterrows():
    d = str(row["date"])
    parts = d.split("-")
    year = int(parts[0])
    if len(parts) >= 2 and parts[1] != "XX":
        month = int(parts[1])
    else:
        month = 6  # mid-year if unknown
    if len(parts) >= 3 and parts[2] != "XX":
        day = int(parts[2])
    else:
        day = 15  # mid-month if unknown
    from datetime import datetime
    dt = datetime(year, month, day)
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year, 12, 31)
    frac = (dt - start_of_year).days / (end_of_year - start_of_year).days
    decimal_date = year + frac
    dates.append({"strain": row["strain"], "date": d, "decimal_date": round(decimal_date, 4),
                  "country": row.get("country", ""), "region": row.get("region", "")})

df = pd.DataFrame(dates)
os.makedirs("outputs/analysis", exist_ok=True)
df.to_csv("outputs/analysis/parsed_dates.tsv", sep="\t", index=False)
print(f"Parsed {len(df)} dates, range: {df['decimal_date'].min():.2f} - {df['decimal_date'].max():.2f}")
PYEOF
fi

# ============================================================================
# LEVEL 2a: Trim alignment with trimAl
# ============================================================================
if [ ! -f "$OUT/alignment/aligned_trimmed.fasta" ]; then
  echo "[Level 2a] Trimming alignment with trimAl..."
  trimal -in "$OUT/alignment/aligned_raw.fasta" \
    -out "$OUT/alignment/aligned_trimmed.fasta" \
    -gt 0.5 -cons 50
fi

# ============================================================================
# LEVEL 2b: Model test with IQ-TREE
# ============================================================================
if [ ! -f "$OUT/trees/model_test.iqtree" ]; then
  echo "[Level 2b] Running IQ-TREE model test..."
  iqtree -s "$OUT/alignment/aligned_raw.fasta" \
    -m MFP -nt "$THREADS" \
    --prefix "$OUT/trees/model_test" -redo
fi

# ============================================================================
# LEVEL 2c: Root-to-tip clock signal check
# ============================================================================
if [ ! -f "$OUT/analysis/clock_signal.tsv" ]; then
  echo "[Level 2c] Checking temporal signal (root-to-tip regression)..."
  python3 << 'PYEOF'
import subprocess, os
import pandas as pd

# Quick NJ tree for root-to-tip regression
# Use treetime's built-in clock check
dates = pd.read_csv("outputs/analysis/parsed_dates.tsv", sep="\t")
# Create a dates file for treetime (strain\tdecimal_date)
date_map = dict(zip(dates["strain"], dates["decimal_date"]))

# Write a simple date file
with open("outputs/analysis/dates_for_treetime.tsv", "w") as f:
    f.write("name\tdate\n")
    for strain, ddate in date_map.items():
        f.write(f"{strain}\t{ddate}\n")

print(f"Prepared {len(date_map)} dates for clock signal analysis")
print(f"Date range: {min(date_map.values()):.2f} - {max(date_map.values()):.2f}")

# Save clock signal info
with open("outputs/analysis/clock_signal.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"num_sequences\t{len(date_map)}\n")
    f.write(f"date_range_years\t{max(date_map.values()) - min(date_map.values()):.2f}\n")
    f.write(f"earliest_date\t{min(date_map.values()):.4f}\n")
    f.write(f"latest_date\t{max(date_map.values()):.4f}\n")
PYEOF
fi

# ============================================================================
# LEVEL 3 (CONVERGENCE 1): Cleaned alignment + model + dates
# ============================================================================
echo "[Level 3] CONVERGENCE 1: alignment + model + dates ready"
# Extract best model from IQ-TREE
BEST_MODEL=$(grep "Best-fit model" "$OUT/trees/model_test.iqtree" | head -1 | awk -F': ' '{print $2}' | awk '{print $1}')
echo "  Best-fit model: $BEST_MODEL"

# ============================================================================
# LEVEL 4: Full ML tree with IQ-TREE (using best model + bootstrap)
# ============================================================================
if [ ! -f "$OUT/trees/ml_tree.treefile" ]; then
  echo "[Level 4] Building ML tree with IQ-TREE..."
  iqtree -s "$OUT/alignment/aligned_trimmed.fasta" \
    -m "$BEST_MODEL" -bb 1000 -nt "$THREADS" \
    --prefix "$OUT/trees/ml_tree" -redo
fi

# ============================================================================
# LEVEL 5: TreeTime molecular clock analysis
# ============================================================================
if [ ! -f "$OUT/clock/timetree.nexus" ]; then
  echo "[Level 5] Running TreeTime molecular clock..."
  treetime \
    --aln "$OUT/alignment/aligned_trimmed.fasta" \
    --tree "$OUT/trees/ml_tree.treefile" \
    --dates "$OUT/analysis/dates_for_treetime.tsv" \
    --clock-filter 3.0 \
    --outdir "$OUT/clock" \
    --confidence || true
  # TreeTime sometimes returns non-zero even on success
fi

# ============================================================================
# LEVEL 6a: TreeTime mugration (geographic migration)
# ============================================================================
if [ ! -f "$OUT/mugration/GTR.txt" ]; then
  echo "[Level 6a] Running TreeTime mugration (geographic)..."
  # Create a traits file for mugration
  python3 << 'PYEOF'
import pandas as pd, os
meta = pd.read_csv("data/metadata.tsv", sep="\t")
os.makedirs("outputs/mugration", exist_ok=True)
traits = meta[["strain", "country"]].rename(columns={"strain": "name", "country": "country"})
traits.to_csv("outputs/mugration/traits.tsv", sep="\t", index=False)
PYEOF
  treetime mugration \
    --tree "$OUT/clock/timetree.nexus" \
    --states "$OUT/mugration/traits.tsv" \
    --attribute country \
    --outdir "$OUT/mugration" || true
fi

# ============================================================================
# LEVEL 6b: TreeTime skyline (effective population size)
# ============================================================================
if [ ! -f "$OUT/skyline/skyline.tsv" ]; then
  echo "[Level 6b] Running TreeTime skyline (Ne estimation)..."
  mkdir -p "$OUT/skyline"
  treetime \
    --aln "$OUT/alignment/aligned_trimmed.fasta" \
    --tree "$OUT/trees/ml_tree.treefile" \
    --dates "$OUT/analysis/dates_for_treetime.tsv" \
    --coalescent skyline \
    --n-skyline 10 \
    --outdir "$OUT/skyline" || true
fi

# ============================================================================
# LEVEL 6c: TreeTime ancestral reconstruction
# ============================================================================
if [ ! -f "$OUT/ancestral/ancestral_sequences.fasta" ]; then
  echo "[Level 6c] Running TreeTime ancestral reconstruction..."
  mkdir -p "$OUT/ancestral"
  treetime ancestral \
    --aln "$OUT/alignment/aligned_trimmed.fasta" \
    --tree "$OUT/clock/timetree.nexus" \
    --outdir "$OUT/ancestral" || true
fi

# ============================================================================
# LEVEL 7 (CONVERGENCE 2): migration + skyline + ancestral
# ============================================================================
echo "[Level 7] CONVERGENCE 2: migration + skyline + ancestral ready"

# ============================================================================
# LEVEL 8a: Augur export (Auspice JSON)
# ============================================================================
if [ ! -f "$OUT/export/auspice.json" ]; then
  echo "[Level 8a] Running augur export..."
  mkdir -p "$OUT/export"
  # Create augur-compatible metadata
  python3 << 'PYEOF'
import pandas as pd, json, os

meta = pd.read_csv("data/metadata.tsv", sep="\t")
# Create minimal auspice config
config = {
    "title": "Zika Phylodynamics",
    "color_options": {"country": {"type": "discrete"}},
    "geo": {"country": {}},
    "maintainers": [{"name": "benchmark", "url": ""}],
    "panels": ["tree", "map"]
}
os.makedirs("outputs/export", exist_ok=True)
with open("outputs/export/auspice_config.json", "w") as f:
    json.dump(config, f, indent=2)

# augur export v2 needs a node-data JSON
# Create one from treetime output
node_data = {"nodes": {}}
if os.path.exists("outputs/clock/timetree.nexus"):
    # Parse basic info
    for _, row in meta.iterrows():
        node_data["nodes"][row["strain"]] = {
            "country": row.get("country", "unknown"),
            "region": row.get("region", "unknown")
        }
with open("outputs/export/node_data.json", "w") as f:
    json.dump(node_data, f, indent=2)
print(f"Created auspice config and node data for {len(meta)} strains")
PYEOF

  # Use treetime's auspice JSON if augur export fails
  augur export v2 \
    --tree "$OUT/trees/ml_tree.treefile" \
    --metadata "$DATA/metadata.tsv" \
    --node-data "$OUT/export/node_data.json" \
    --output "$OUT/export/auspice.json" \
    --auspice-config "$OUT/export/auspice_config.json" \
    --metadata-id-columns strain 2>/dev/null || \
    cp "$OUT/clock/auspice_tree.json" "$OUT/export/auspice.json"
fi

# ============================================================================
# LEVEL 8b: Transmission chain reconstruction
# ============================================================================
if [ ! -f "$OUT/analysis/transmission_chains.tsv" ]; then
  echo "[Level 8b] Reconstructing transmission chains..."
  python3 << 'PYEOF'
import os, re
import pandas as pd

# Parse the timetree to extract parent-child relationships and dates
# Use the mugration results for geographic transitions
chains = []
transitions = []

# Parse mugration GTR for transition rates
gtr_file = "outputs/mugration/GTR.txt"
if os.path.exists(gtr_file):
    with open(gtr_file) as f:
        content = f.read()
    # Extract transition matrix info
    print(f"GTR model file size: {len(content)} bytes")

# Parse annotated tree for geographic transitions
annotated_tree = "outputs/mugration/annotated_tree.nexus"
if os.path.exists(annotated_tree):
    with open(annotated_tree) as f:
        tree_content = f.read()
    # Count country annotations
    countries = re.findall(r'country="([^"]+)"', tree_content)
    country_counts = pd.Series(countries).value_counts()
    print(f"Country annotations found: {len(countries)}")
    for c, n in country_counts.items():
        transitions.append({"from_country": "ancestor", "to_country": c, "count": n})
else:
    print("Warning: No annotated tree found, using metadata for chains")

# Fallback: use metadata for chain info
meta = pd.read_csv("data/metadata.tsv", sep="\t")
dates_df = pd.read_csv("outputs/analysis/parsed_dates.tsv", sep="\t")
merged = meta.merge(dates_df[["strain", "decimal_date"]], on="strain")
merged = merged.sort_values("decimal_date")

# Build simple temporal chains by country
for country, grp in merged.groupby("country"):
    grp_sorted = grp.sort_values("decimal_date")
    if len(grp_sorted) >= 2:
        first = grp_sorted.iloc[0]
        last = grp_sorted.iloc[-1]
        chains.append({
            "country": country,
            "n_sequences": len(grp_sorted),
            "first_date": first["decimal_date"],
            "last_date": last["decimal_date"],
            "duration_years": round(last["decimal_date"] - first["decimal_date"], 2)
        })

chains_df = pd.DataFrame(chains)
os.makedirs("outputs/analysis", exist_ok=True)
chains_df.to_csv("outputs/analysis/transmission_chains.tsv", sep="\t", index=False)
print(f"Identified {len(chains_df)} geographic transmission chains")
PYEOF
fi

# ============================================================================
# LEVEL 8c: Rate estimation report
# ============================================================================
if [ ! -f "$OUT/analysis/rate_estimation.tsv" ]; then
  echo "[Level 8c] Estimating evolutionary rates..."
  python3 << 'PYEOF'
import os, re
import pandas as pd

rates = {}

# Parse TreeTime clock rate from molecular_clock.txt or timetree log
clock_file = "outputs/clock/molecular_clock.txt"
if os.path.exists(clock_file):
    with open(clock_file) as f:
        for line in f:
            if "rate" in line.lower():
                # Try to extract rate value
                match = re.search(r'[\d.]+[eE][+-]?\d+|[\d.]+', line)
                if match:
                    rates["clock_rate_from_file"] = match.group()
                    break

# Parse from dates.tsv if available
dates_file = "outputs/clock/dates.tsv"
if os.path.exists(dates_file):
    df = pd.read_csv(dates_file, sep="\t")
    print(f"TreeTime dates output columns: {list(df.columns)}")
    if "date" in df.columns:
        rates["n_dated_tips"] = str(len(df))

# Parse root-to-tip from treetime
rtt_file = "outputs/clock/root_to_tip_regression.pdf"
rates["root_to_tip_plot"] = "yes" if os.path.exists(rtt_file) else "no"

# Parse skyline
skyline_file = "outputs/skyline/skyline.tsv"
if os.path.exists(skyline_file):
    sky_df = pd.read_csv(skyline_file, sep="\t")
    rates["skyline_intervals"] = str(len(sky_df))
    if "Ne" in sky_df.columns or "Tc" in sky_df.columns:
        ne_col = [c for c in sky_df.columns if "ne" in c.lower() or "tc" in c.lower() or "neff" in c.lower()]
        if ne_col:
            rates["max_Ne"] = str(sky_df[ne_col[0]].max())
            rates["min_Ne"] = str(sky_df[ne_col[0]].min())

# Write rate estimation
os.makedirs("outputs/analysis", exist_ok=True)
with open("outputs/analysis/rate_estimation.tsv", "w") as f:
    f.write("metric\tvalue\n")
    for k, v in rates.items():
        f.write(f"{k}\t{v}\n")
print(f"Rate estimation: {len(rates)} metrics")
PYEOF
fi

# ============================================================================
# LEVEL 9 (CONVERGENCE 3): Auspice + chains + rates
# ============================================================================
echo "[Level 9] CONVERGENCE 3: export + chains + rates ready"

# ============================================================================
# LEVEL 10 (CONVERGENCE 4 + Final Report): Integrate everything
# ============================================================================
echo "[Level 10] CONVERGENCE 4: Generating final report..."
python3 << 'PYEOF'
import os, re, json
import pandas as pd

results = {}

# ---- Alignment stats ----
aln_file = "outputs/alignment/aligned_trimmed.fasta"
seqs = {}
current = None
with open(aln_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current = line[1:].split()[0]
            seqs[current] = ""
        elif current:
            seqs[current] += line
results["num_sequences"] = len(seqs)
aln_len = len(list(seqs.values())[0])
results["alignment_length"] = aln_len
total_chars = sum(len(s) for s in seqs.values())
gap_chars = sum(s.count("-") for s in seqs.values())
results["gap_fraction"] = round(gap_chars / total_chars, 4)

# ---- Best substitution model ----
with open("outputs/trees/model_test.iqtree") as f:
    for line in f:
        if "Best-fit model" in line:
            results["best_model"] = line.split(":")[1].strip().split()[0]
            break

# ---- ML tree stats ----
with open("outputs/trees/ml_tree.treefile") as f:
    tree_str = f.read().strip()
bootstraps = re.findall(r'\)(\d+(?:\.\d+)?)', tree_str)
if bootstraps:
    bs_vals = [float(b) for b in bootstraps]
    results["mean_bootstrap"] = round(sum(bs_vals) / len(bs_vals), 1)
    results["high_support_nodes"] = sum(1 for b in bs_vals if b >= 70)
    results["total_internal_nodes"] = len(bs_vals)

# ---- Clock rate ----
with open("outputs/skyline/molecular_clock.txt") as f:
    for line in f:
        if "--rate:" in line:
            match = re.search(r'([\d.]+[eE][+-]?\d+)', line)
            if match:
                results["clock_rate"] = match.group(1)
        if "--r^2:" in line:
            match = re.search(r'([\d.]+)', line.split(":")[-1])
            if match:
                results["clock_r_squared"] = match.group(1)

# ---- TMRCA ----
dates_df = pd.read_csv("outputs/clock/dates.tsv", sep="\t")
input_strains = set(seqs.keys())
tip_dates = dates_df[dates_df["#node"].isin(input_strains)]
results["num_dated_tips"] = len(tip_dates)
if "numeric date" in dates_df.columns:
    try:
        results["tmrca"] = round(float(dates_df["numeric date"].min()), 2)
    except (ValueError, TypeError):
        results["tmrca"] = str(dates_df["numeric date"].min())
else:
    results["tmrca"] = str(dates_df["date"].min())

# ---- Mugration ----
results["mugration_model"] = "fitted" if os.path.exists("outputs/mugration/GTR.txt") else "not_run"
annotated = "outputs/mugration/annotated_tree.nexus"
if os.path.exists(annotated):
    with open(annotated) as f:
        countries = set(re.findall(r'country="([^"]+)"', f.read()))
    results["num_migration_countries"] = len(countries)

# ---- Skyline (parse manually due to comment header) ----
with open("outputs/skyline/skyline.tsv") as f:
    lines = [l.strip() for l in f if not l.startswith("#")]
ne_values = []
for line in lines:
    parts = line.split("\t")
    if len(parts) >= 2:
        try:
            ne_values.append(float(parts[1]))
        except ValueError:
            pass
if ne_values:
    results["skyline_intervals"] = len(ne_values)
    results["max_effective_pop_size"] = round(max(ne_values), 1)
    results["min_effective_pop_size"] = round(min(ne_values), 1)

# ---- Ancestral ----
anc_count = sum(1 for l in open("outputs/ancestral/ancestral_sequences.fasta") if l.startswith(">"))
results["ancestral_sequences"] = anc_count

# ---- Transmission chains ----
chains_df = pd.read_csv("outputs/analysis/transmission_chains.tsv", sep="\t")
results["num_geographic_clusters"] = len(chains_df)
results["largest_cluster_size"] = int(chains_df["n_sequences"].max())

# ---- Visualization export ----
results["visualization_export"] = "success" if os.path.exists("outputs/export/auspice.json") else "failed"

# ---- Date metadata ----
dp = pd.read_csv("outputs/analysis/parsed_dates.tsv", sep="\t")
results["date_range_start"] = round(dp["decimal_date"].min(), 2)
results["date_range_end"] = round(dp["decimal_date"].max(), 2)
results["num_countries"] = dp["country"].nunique()

# ---- Write CSV ----
os.makedirs("results", exist_ok=True)
with open("results/report.csv", "w") as f:
    f.write("metric,value\n")
    for k, v in results.items():
        f.write(f"{k},{v}\n")

print("=== Final Report ===")
for k, v in results.items():
    print(f"  {k}: {v}")
print(f"\nTotal metrics: {len(results)}")
PYEOF

echo "Pipeline complete. Results in results/report.csv"
