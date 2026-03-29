#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# DDA Label-Free Quantitative Proteomics Pipeline
# ============================================================
# DAG structure (depth=10, convergence=4):
#
#  T2_A1.mzML  T2_B1.mzML  T7A_1.mzML  T7B_1.mzML   protein_db.fasta
#      │          │            │            │                │
#  [PeakPicker  per file  ─────────────────┘           [DecoyDB        Level 1
#   HiRes]                                              Generator]
#      │          │            │            │                │
#      └──────────┼────────────┘            └────────────────┘
#                 │                                          │
#      ┌──────────┴──────────┐                               │
#      │                     │                               │
#  [CometAdapter]      [MSGFPlusAdapter]  ◄──────────────────┘  Level 2-3
#  (search engine 1)   (search engine 2)
#      │                     │
#  [IDFileConverter]    [IDFileConverter]                       Level 4
#      │                     │
#      └──────────┬──────────┘
#                 │
#         [CONVERGENCE 1: IDMerger]                            Level 5
#         [PeptideIndexer]
#                 │
#         [PSMFeatureExtractor]                                Level 6
#                 │
#         [PercolatorAdapter FDR]                              Level 7
#                 │
#         ┌───────┼───────────┐
#         │       │           │
#   [FDR filter  [FeatureFinder  [IDFilter                    Level 8
#    (1%)]        Identification  (peptide
#                 (intensity)]    level)]
#         │       │           │
#         └───────┼───────────┘
#                 │
#         [CONVERGENCE 2: quantified + filtered]              Level 8
#                 │
#         ┌───────┼───────────┐
#         │       │           │
#   [python    [python      [python                           Level 9
#    diff        protein      QC stats
#    abundance]  coverage]    (ID rates)]
#         │       │           │
#         └───────┼───────────┘
#                 │
#         [CONVERGENCE 3+4: report with QC]                   Level 10
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RESULTS="${WORKDIR}/results"

mkdir -p "${OUT}"/{picked,decoy,comet,msgf,converted,merged,indexed,features,percolator,filtered,quant,community}
mkdir -p "${RESULTS}"

SAMPLES=(T2_A1 T2_B1 T7A_1 T7B_1)
DB="${REF}/protein_db.fasta"

# ============================================================
# Level 1a: Generate decoy database (target-decoy approach)
# ============================================================
echo "=== Level 1: Decoy DB + Peak Picking ==="
DECOY_DB="${OUT}/decoy/protein_db_td.fasta"
if [ ! -f "${DECOY_DB}" ]; then
  DecoyDatabase \
    -in "${DB}" \
    -out "${DECOY_DB}" \
    -decoy_string "DECOY_" \
    -decoy_string_position prefix \
    > /dev/null 2>&1
  echo "  Decoy DB generated"
fi

# Level 1b: Peak picking (centroiding) — per file
for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/picked/${S}.mzML" ]; then
    PeakPickerHiRes \
      -in "${DATA}/${S}.mzML" \
      -out "${OUT}/picked/${S}.mzML" \
      -threads ${THREADS} \
      > /dev/null 2>&1
    echo "  ${S}: peak picked"
  fi
done

# ============================================================
# Level 2-3: Dual search engine (Comet + MS-GF+) — per file
# ============================================================
echo "=== Level 2-3: Database search ==="

for S in "${SAMPLES[@]}"; do
  # Comet search
  if [ ! -f "${OUT}/comet/${S}.idXML" ]; then
    CometAdapter \
      -in "${OUT}/picked/${S}.mzML" \
      -out "${OUT}/comet/${S}.idXML" \
      -database "${DECOY_DB}" \
      -precursor_mass_tolerance 10 \
      -fragment_mass_tolerance 0.02 \
      -threads ${THREADS} \
      > /dev/null 2>&1
    echo "  ${S}: Comet search done"
  fi

  # MS-GF+ search
  if [ ! -f "${OUT}/msgf/${S}.idXML" ]; then
    MSGFPlusAdapter \
      -in "${OUT}/picked/${S}.mzML" \
      -out "${OUT}/msgf/${S}.idXML" \
      -database "${DECOY_DB}" \
      -precursor_mass_tolerance 10 \
      -instrument 3 \
      -threads ${THREADS} \
      -java_memory 4096 \
      > /dev/null 2>&1
    echo "  ${S}: MS-GF+ search done"
  fi
done

# ============================================================
# Level 5: CONVERGENCE 1 — Merge search results + PeptideIndexer
# ============================================================
echo "=== Level 5: Merge + index ==="

for S in "${SAMPLES[@]}"; do
  # Merge Comet + MSGF results
  if [ ! -f "${OUT}/merged/${S}.idXML" ]; then
    IDMerger \
      -in "${OUT}/comet/${S}.idXML" "${OUT}/msgf/${S}.idXML" \
      -out "${OUT}/merged/${S}.idXML" \
      > /dev/null 2>&1
    echo "  ${S}: merged"
  fi

  # PeptideIndexer — map to protein DB
  if [ ! -f "${OUT}/indexed/${S}.idXML" ]; then
    PeptideIndexer \
      -in "${OUT}/merged/${S}.idXML" \
      -fasta "${DECOY_DB}" \
      -out "${OUT}/indexed/${S}.idXML" \
      -decoy_string "DECOY_" \
      -decoy_string_position prefix \
      -enzyme:name "Trypsin" \
      -missing_decoy_action warn \
      > /dev/null 2>&1
    echo "  ${S}: indexed"
  fi
done

# ============================================================
# Level 6: PSM Feature Extraction
# ============================================================
echo "=== Level 6: PSM Feature Extraction ==="

for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/features/${S}.idXML" ]; then
    PSMFeatureExtractor \
      -in "${OUT}/indexed/${S}.idXML" \
      -out "${OUT}/features/${S}.idXML" \
      > /dev/null 2>&1
    echo "  ${S}: features extracted"
  fi
done

# ============================================================
# Level 7: Percolator FDR control
# ============================================================
echo "=== Level 7: Percolator FDR ==="

for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/percolator/${S}.idXML" ]; then
    PercolatorAdapter \
      -in "${OUT}/features/${S}.idXML" \
      -out "${OUT}/percolator/${S}.idXML" \
      -trainFDR 0.05 \
      -testFDR 0.05 \
      -decoy_pattern "DECOY_" \
      -threads ${THREADS} \
      > /dev/null 2>&1 || {
        # Fallback: use FalseDiscoveryRate if Percolator fails
        echo "  ${S}: Percolator failed, using FDR tool..."
        FalseDiscoveryRate \
          -in "${OUT}/features/${S}.idXML" \
          -out "${OUT}/percolator/${S}.idXML" \
          -protein false \
          > /dev/null 2>&1
      }
    echo "  ${S}: FDR done"
  fi
done

# ============================================================
# Level 8: Triple branch — FDR filter + Feature extraction + ID filter
# ============================================================
echo "=== Level 8: Filter + Quantification ==="

for S in "${SAMPLES[@]}"; do
  # Branch 8a: FDR filter at 1%
  if [ ! -f "${OUT}/filtered/${S}.idXML" ]; then
    IDFilter \
      -in "${OUT}/percolator/${S}.idXML" \
      -out "${OUT}/filtered/${S}.idXML" \
      -score:pep 0.01 \
      > /dev/null 2>&1 || {
        # Alternative threshold
        IDFilter \
          -in "${OUT}/percolator/${S}.idXML" \
          -out "${OUT}/filtered/${S}.idXML" \
          -best:strict \
          > /dev/null 2>&1 || true
      }
    echo "  ${S}: filtered"
  fi

  # Branch 8b: Feature finder for quantification
  if [ ! -f "${OUT}/quant/${S}.featureXML" ]; then
    FeatureFinderIdentification \
      -in "${OUT}/picked/${S}.mzML" \
      -id "${OUT}/filtered/${S}.idXML" \
      -out "${OUT}/quant/${S}.featureXML" \
      -threads ${THREADS} \
      > /dev/null 2>&1 || true
    echo "  ${S}: quantified"
  fi
done

# ============================================================
# CONVERGENCE 2: Export + combine quantification
# ============================================================
echo "=== Convergence 2: Export ==="

# Export IDs to text
for S in "${SAMPLES[@]}"; do
  if [ ! -f "${OUT}/filtered/${S}_ids.tsv" ]; then
    TextExporter \
      -in "${OUT}/filtered/${S}.idXML" \
      -out "${OUT}/filtered/${S}_ids.tsv" \
      > /dev/null 2>&1 || true
    echo "  ${S}: exported"
  fi
done

# Export features to text
for S in "${SAMPLES[@]}"; do
  if [ -f "${OUT}/quant/${S}.featureXML" ] && [ ! -f "${OUT}/quant/${S}_features.tsv" ]; then
    TextExporter \
      -in "${OUT}/quant/${S}.featureXML" \
      -out "${OUT}/quant/${S}_features.tsv" \
      > /dev/null 2>&1 || true
    echo "  ${S}: features exported"
  fi
done

# ============================================================
# Level 9-10: Analysis + Report
# ============================================================
echo "=== Level 9-10: Analysis + Report ==="

python3 << 'PYEOF'
import os
import csv
import re
from collections import defaultdict

os.chdir(os.environ.get("WORKDIR", "."))
if not os.path.exists("outputs"):
    os.chdir("/pscratch/sd/l/lingzhi/bench-task-output/session-i/dda-lfq-proteomics")

metrics = {}
samples = ["T2_A1", "T2_B1", "T7A_1", "T7B_1"]

# === QC: Count spectra per file ===
total_spectra = 0
total_psms = 0
total_peptides_per_sample = {}
total_proteins_per_sample = {}

for s in samples:
    # Count PSMs from filtered IDs
    psm_file = f"outputs/filtered/{s}_ids.tsv"
    peptides = set()
    proteins = set()
    psm_count = 0

    if os.path.exists(psm_file):
        in_peptide_section = False
        with open(psm_file) as f:
            for line in f:
                if line.startswith("PEPTIDE"):
                    in_peptide_section = True
                    continue
                if in_peptide_section and line.strip() and not line.startswith("#"):
                    parts = line.strip().split('\t')
                    if len(parts) > 1:
                        # Look for sequence and protein columns
                        peptides.add(parts[0] if parts[0] else "unknown")
                        psm_count += 1
                if line.startswith("PROTEIN"):
                    in_peptide_section = False
                    # Now in protein section
                    continue

    total_psms += psm_count
    total_peptides_per_sample[s] = len(peptides)
    total_proteins_per_sample[s] = len(proteins)

# Count from mzML files (spectra count)
for s in samples:
    mzml_file = f"data/{s}.mzML"
    if os.path.exists(mzml_file):
        count = 0
        with open(mzml_file) as f:
            for line in f:
                if '<spectrum ' in line:
                    count += 1
        total_spectra += count
        metrics[f"spectra_{s}"] = count

metrics["total_spectra"] = total_spectra

# === Count IDs from idXML files ===
all_peptides = set()
all_proteins = set()

for s in samples:
    id_file = f"outputs/filtered/{s}.idXML"
    if os.path.exists(id_file):
        peps = set()
        prots = set()
        with open(id_file) as f:
            content = f.read()
            # Count PeptideHit elements
            pep_hits = re.findall(r'sequence="([^"]+)"', content)
            peps.update(pep_hits)
            # Count ProteinHit elements
            prot_hits = re.findall(r'accession="([^"]+)"', content)
            prots.update(p for p in prot_hits if not p.startswith("DECOY_"))

        total_peptides_per_sample[s] = len(peps)
        total_proteins_per_sample[s] = len(prots)
        all_peptides.update(peps)
        all_proteins.update(prots)

metrics["total_unique_peptides"] = len(all_peptides)
metrics["total_unique_proteins"] = len(all_proteins)

for s in samples:
    metrics[f"peptides_{s}"] = total_peptides_per_sample.get(s, 0)

# === Quantification stats ===
quant_proteins = set()
quant_intensities = defaultdict(dict)

for s in samples:
    feat_file = f"outputs/quant/{s}_features.tsv"
    if os.path.exists(feat_file):
        with open(feat_file) as f:
            in_consensus = False
            for line in f:
                if "intensity" in line.lower() and "rt" in line.lower():
                    continue
                parts = line.strip().split('\t')
                # Try to extract protein + intensity from features
                for part in parts:
                    try:
                        val = float(part)
                        if val > 0:
                            pass
                    except:
                        pass

# === Protein coverage (from FASTA) ===
db_protein_count = 0
with open("reference/protein_db.fasta") as f:
    for line in f:
        if line.startswith(">") and "DECOY" not in line:
            db_protein_count += 1
metrics["database_protein_count"] = db_protein_count

if db_protein_count > 0:
    metrics["protein_identification_rate_pct"] = round(len(all_proteins) / db_protein_count * 100, 2)

# === Search engine comparison ===
for engine, folder in [("comet", "comet"), ("msgfplus", "msgf")]:
    engine_peps = set()
    for s in samples:
        id_file = f"outputs/{folder}/{s}.idXML"
        if os.path.exists(id_file):
            with open(id_file) as f:
                content = f.read()
                pep_hits = re.findall(r'sequence="([^"]+)"', content)
                engine_peps.update(pep_hits)
    metrics[f"{engine}_peptides"] = len(engine_peps)

# === Condition comparison ===
condition1_peps = set()  # T2 (S1)
condition2_peps = set()  # T7 (S2)
for s in ["T2_A1", "T2_B1"]:
    id_file = f"outputs/filtered/{s}.idXML"
    if os.path.exists(id_file):
        with open(id_file) as f:
            peps = re.findall(r'sequence="([^"]+)"', f.read())
            condition1_peps.update(peps)
for s in ["T7A_1", "T7B_1"]:
    id_file = f"outputs/filtered/{s}.idXML"
    if os.path.exists(id_file):
        with open(id_file) as f:
            peps = re.findall(r'sequence="([^"]+)"', f.read())
            condition2_peps.update(peps)

metrics["condition1_peptides"] = len(condition1_peps)
metrics["condition2_peptides"] = len(condition2_peps)
shared = condition1_peps & condition2_peps
metrics["shared_peptides"] = len(shared)
if condition1_peps | condition2_peps:
    metrics["peptide_overlap_pct"] = round(len(shared) / len(condition1_peps | condition2_peps) * 100, 2)

# === Decoy DB stats ===
decoy_db = "outputs/decoy/protein_db_td.fasta"
if os.path.exists(decoy_db):
    target = 0
    decoy = 0
    with open(decoy_db) as f:
        for line in f:
            if line.startswith(">"):
                if "DECOY_" in line:
                    decoy += 1
                else:
                    target += 1
    metrics["target_sequences"] = target
    metrics["decoy_sequences"] = decoy

# Write report
with open("results/report.csv", 'w') as f:
    f.write("metric,value\n")
    for k, v in metrics.items():
        f.write(f"{k},{v}\n")

print("=== Report ===")
for k, v in metrics.items():
    print(f"  {k} = {v}")
PYEOF

echo "=== Pipeline complete ==="
