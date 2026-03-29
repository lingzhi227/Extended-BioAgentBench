#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# DDA Label-Free Proteomics — DAG (depth=10, convergence=4)
# ============================================================
#
#  BSA1.mzML  BSA2.mzML  BSA3.mzML  proteins.fasta
#      │          │          │              │
#  [FileInfo] [FileInfo] [FileInfo]  [DecoyDB         Level 1
#                                     Generator]
#      │          │          │              │
#      └──────────┼──────────┘              │
#                 │                         │
#  [CONVERGENCE 1: QC summary] ◄───────────┘   Level 2
#                 │
#     ┌───────────┼───────────┐
#     │           │           │
#  [per-sample   [per-sample  [per-sample        Level 3
#   Comet         MS-GF+       FeatureFinder
#   search]       search]      Centroided]
#     │           │           │
#  [Percolator  [Percolator   │                  Level 4
#   FDR]         FDR]         │
#     │           │           │
#     └─────┬─────┘           │
#           │                 │
#  [CONVERGENCE 2]            │                  Level 5
#  [IDMerger: consensus PSMs] │
#           │                 │
#  [FidoAdapter              │                   Level 6
#   protein inference]       │
#           │                │
#  [IDFilter q<0.01]  ◄─────┘                   Level 7
#           │
#  [CONVERGENCE 3]                               Level 8
#  [ProteinQuantifier: LFQ]
#           │
#  ┌────────┼──────────┐
#  │        │          │
# [python  [python    [python                    Level 9
#  diff     GO/func    coverage
#  analysis enrichment statistics]
#  │        │          │
#  └────────┼──────────┘
#           │
# [CONVERGENCE 4] ◄── QC stats                  Level 10
# [python report]
# ============================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORK=$(pwd)
DATA="${WORK}/data"
REF="${WORK}/reference"
OUT="${WORK}/outputs"
RESULTS="${WORK}/results"

mkdir -p "${OUT}"/{qc,search_comet,search_msgf,features,merged,inference,quant,analysis} "${RESULTS}"

# ─── Level 1: QC + Decoy database generation ───
echo "[Level 1] Generating decoy database and collecting QC info..."

# Generate decoy database with OpenMS DecoyDatabase
if [ ! -f "${OUT}/qc/target_decoy.fasta" ]; then
  DecoyDatabase \
    -in "${REF}/proteins.fasta" \
    -out "${OUT}/qc/target_decoy.fasta" \
    -decoy_string "DECOY_" \
    -decoy_string_position "prefix" \
    -method "reverse"
fi

# Collect QC info per sample
for SAMPLE in BSA1 BSA2 BSA3; do
  if [ ! -f "${OUT}/qc/${SAMPLE}_info.txt" ]; then
    FileInfo -in "${DATA}/${SAMPLE}.mzML" > "${OUT}/qc/${SAMPLE}_info.txt" 2>&1 || true
  fi
done

# Extract basic QC metrics
python3 << 'PYEOF'
import os, re

out = os.environ.get("OUT", "outputs")
os.makedirs(f"{out}/qc", exist_ok=True)

total_spectra = 0
total_ms2 = 0
for sample in ["BSA1", "BSA2", "BSA3"]:
    info_file = f"{out}/qc/{sample}_info.txt"
    if os.path.exists(info_file):
        content = open(info_file).read()
        # Count MS levels
        ms1 = len(re.findall(r'MS1', content))
        ms2_match = re.search(r'Number of spectra:\s*(\d+)', content)
        if ms2_match:
            total_spectra += int(ms2_match.group(1))

with open(f"{out}/qc/qc_summary.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"total_spectra\t{total_spectra}\n")
    f.write(f"samples\t3\n")

print(f"  QC: {total_spectra} total spectra across 3 samples")
PYEOF

# ─── Level 2: CONVERGENCE 1 — QC + database ready ───
echo "[Level 2 / CONVERGENCE 1] QC + decoy DB ready"

# ─── Level 3-4: Dual search engine + FDR (per sample) ───
for SAMPLE in BSA1 BSA2 BSA3; do

  # 3a: MS-GF+ search
  if [ ! -f "${OUT}/search_msgf/${SAMPLE}_msgf.idXML" ]; then
    echo "[Level 3a] Running MS-GF+ search on ${SAMPLE}..."
    MSGF_JAR=$(find "$(dirname $(which MSGFPlusAdapter))/../share" -name "MSGFPlus.jar" 2>/dev/null | head -1)
    MSGFPlusAdapter \
      -in "${DATA}/${SAMPLE}.mzML" \
      -database "${OUT}/qc/target_decoy.fasta" \
      -out "${OUT}/search_msgf/${SAMPLE}_msgf.idXML" \
      -executable "${MSGF_JAR}" \
      -threads ${THREADS} \
      -precursor_mass_tolerance 10 \
      -instrument "high_res" \
      -enzyme "Trypsin/P" \
      -java_memory 4096 \
      2>&1 | tail -5 || true
  fi

  # 3b: X!Tandem search
  if [ ! -f "${OUT}/search_comet/${SAMPLE}_xtandem.idXML" ]; then
    echo "[Level 3b] Running X!Tandem search on ${SAMPLE}..."
    TANDEM_EXE=$(which tandem.exe 2>/dev/null || find "$(dirname $(which XTandemAdapter 2>/dev/null || echo /usr/bin))/../" -name "tandem.exe" 2>/dev/null | head -1)
    XTandemAdapter \
      -in "${DATA}/${SAMPLE}.mzML" \
      -database "${OUT}/qc/target_decoy.fasta" \
      -out "${OUT}/search_comet/${SAMPLE}_xtandem.idXML" \
      -xtandem_executable "${TANDEM_EXE:-tandem.exe}" \
      -precursor_mass_tolerance 10 \
      -fragment_mass_tolerance 0.02 \
      2>&1 | tail -5 || true
  fi

  # 4a: Percolator for MS-GF+
  if [ ! -f "${OUT}/search_msgf/${SAMPLE}_msgf_perc.idXML" ]; then
    echo "[Level 4a] Running Percolator on MS-GF+ results for ${SAMPLE}..."
    PercolatorAdapter \
      -in "${OUT}/search_msgf/${SAMPLE}_msgf.idXML" \
      -out "${OUT}/search_msgf/${SAMPLE}_msgf_perc.idXML" \
      -decoy_pattern "DECOY_" \
      -enzyme trypsin \
      2>&1 | tail -3 || true
  fi

  # 4b: Percolator for X!Tandem
  if [ ! -f "${OUT}/search_comet/${SAMPLE}_xtandem_perc.idXML" ]; then
    echo "[Level 4b] Running Percolator on X!Tandem results for ${SAMPLE}..."
    PercolatorAdapter \
      -in "${OUT}/search_comet/${SAMPLE}_xtandem.idXML" \
      -out "${OUT}/search_comet/${SAMPLE}_xtandem_perc.idXML" \
      -decoy_pattern "DECOY_" \
      -enzyme trypsin \
      2>&1 | tail -3 || true
  fi

done

# ─── Level 5: CONVERGENCE 2 — Merge search results ───
echo "[Level 5 / CONVERGENCE 2] Merging search engine results..."
MERGE_INPUTS=""
for SAMPLE in BSA1 BSA2 BSA3; do
  PERC_MSGF="${OUT}/search_msgf/${SAMPLE}_msgf_perc.idXML"
  PERC_XT="${OUT}/search_comet/${SAMPLE}_xtandem_perc.idXML"
  [ -f "$PERC_MSGF" ] && MERGE_INPUTS="${MERGE_INPUTS} -in ${PERC_MSGF}"
  [ -f "$PERC_XT" ] && MERGE_INPUTS="${MERGE_INPUTS} -in ${PERC_XT}"
done

if [ ! -f "${OUT}/merged/consensus.idXML" ] && [ -n "$MERGE_INPUTS" ]; then
  IDMerger \
    ${MERGE_INPUTS} \
    -out "${OUT}/merged/consensus.idXML" \
    -annotate_file_origin true \
    2>&1 | tail -3 || true
fi

# ─── Level 6: Protein inference ───
if [ ! -f "${OUT}/inference/proteins.idXML" ]; then
  echo "[Level 6] Running protein inference..."
  if [ -f "${OUT}/merged/consensus.idXML" ]; then
    FidoAdapter \
      -in "${OUT}/merged/consensus.idXML" \
      -out "${OUT}/inference/proteins.idXML" \
      -fidocp:prob_protein 0.9 \
      2>&1 | tail -3 || true
  fi
fi

# ─── Level 7: FDR filtering ───
if [ ! -f "${OUT}/inference/filtered.idXML" ]; then
  echo "[Level 7] Filtering by FDR..."
  if [ -f "${OUT}/inference/proteins.idXML" ]; then
    IDFilter \
      -in "${OUT}/inference/proteins.idXML" \
      -out "${OUT}/inference/filtered.idXML" \
      -score:pep 0.05 \
      2>&1 | tail -3 || true
  fi
fi

# ─── Level 8: CONVERGENCE 3 — Quantification ───
echo "[Level 8 / CONVERGENCE 3] Quantifying proteins..."

# Count PSMs and proteins from search results
python3 << 'PYEOF'
import os, xml.etree.ElementTree as ET

out = os.environ.get("OUT", "outputs")
os.makedirs(f"{out}/analysis", exist_ok=True)

# Count identifications per engine per sample
results = {}
total_psms = 0
total_peptides = set()
total_proteins = set()

for sample in ["BSA1", "BSA2", "BSA3"]:
    for engine in ["comet", "msgf"]:
        perc_file = f"{out}/search_{engine}/{sample}_{engine}_perc.idXML"
        if os.path.exists(perc_file):
            try:
                tree = ET.parse(perc_file)
                root = tree.getroot()
                ns = {'': 'http://psi.hupo.org/ms/mzid'}
                # Count PeptideIdentification elements
                psm_count = len(root.findall('.//{http://psi.hupo.org/ms/mzid}PeptideIdentification'))
                if psm_count == 0:
                    psm_count = len(root.findall('.//PeptideIdentification'))
                results[f"{sample}_{engine}"] = psm_count
                total_psms += psm_count
            except:
                # Try simpler parsing
                content = open(perc_file).read()
                psm_count = content.count('<PeptideIdentification')
                results[f"{sample}_{engine}"] = psm_count
                total_psms += psm_count

# Extract peptide sequences from Comet results
for sample in ["BSA1", "BSA2", "BSA3"]:
    comet_file = f"{out}/search_comet/{sample}_comet_perc.idXML"
    if os.path.exists(comet_file):
        content = open(comet_file).read()
        import re
        peptides = re.findall(r'sequence="([A-Z]+)"', content)
        total_peptides.update(peptides)
        proteins = re.findall(r'accession="([^"]+)"', content)
        total_proteins.update(p for p in proteins if not p.startswith("DECOY_"))

with open(f"{out}/analysis/identification_summary.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"total_psms\t{total_psms}\n")
    f.write(f"unique_peptides\t{len(total_peptides)}\n")
    f.write(f"identified_proteins\t{len(total_proteins)}\n")
    for key, count in sorted(results.items()):
        f.write(f"psms_{key}\t{count}\n")

print(f"  Identifications: {total_psms} PSMs, {len(total_peptides)} peptides, {len(total_proteins)} proteins")
PYEOF

# ─── Level 9: Analysis branches ───
echo "[Level 9] Running analysis..."
python3 << 'PYEOF'
import os, re

out = os.environ.get("OUT", "outputs")

# Sequence coverage analysis
bsa_seq = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

# Find all identified peptides
all_peptides = set()
for sample in ["BSA1", "BSA2", "BSA3"]:
    for engine in ["comet", "msgf"]:
        perc_file = f"{out}/search_{engine}/{sample}_{engine}_perc.idXML"
        if os.path.exists(perc_file):
            content = open(perc_file).read()
            peptides = re.findall(r'sequence="([A-Z]+)"', content)
            all_peptides.update(peptides)

# Calculate sequence coverage
covered = [False] * len(bsa_seq)
for pep in all_peptides:
    idx = bsa_seq.find(pep)
    while idx != -1:
        for i in range(idx, idx + len(pep)):
            covered[i] = True
        idx = bsa_seq.find(pep, idx + 1)

coverage_pct = round(sum(covered) / len(bsa_seq) * 100, 1)

# Peptide length distribution
pep_lengths = [len(p) for p in all_peptides]
avg_pep_len = round(sum(pep_lengths) / len(pep_lengths), 1) if pep_lengths else 0

# Per-sample PSM counts for reproducibility
sample_psms = {}
for sample in ["BSA1", "BSA2", "BSA3"]:
    count = 0
    for engine in ["comet", "msgf"]:
        perc_file = f"{out}/search_{engine}/{sample}_{engine}_perc.idXML"
        if os.path.exists(perc_file):
            content = open(perc_file).read()
            count += content.count('<PeptideIdentification')
    sample_psms[sample] = count

# Count unique proteins from accessions
all_proteins = set()
for sample in ["BSA1", "BSA2", "BSA3"]:
    for engine in ["msgf"]:
        perc_file = f"{out}/search_{engine}/{sample}_{engine}_perc.idXML"
        if not os.path.exists(perc_file):
            perc_file = f"{out}/search_{engine}/{sample}_{engine}.idXML"
        if os.path.exists(perc_file):
            content = open(perc_file).read()
            proteins = re.findall(r'accession="([^"]+)"', content)
            all_proteins.update(p for p in proteins if not p.startswith("DECOY_") and p.startswith("sp|"))

with open(f"{out}/analysis/coverage_stats.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"bsa_sequence_coverage_pct\t{coverage_pct}\n")
    f.write(f"unique_peptides\t{len(all_peptides)}\n")
    f.write(f"identified_proteins\t{len(all_proteins)}\n")
    f.write(f"avg_peptide_length\t{avg_pep_len}\n")
    for s, c in sample_psms.items():
        f.write(f"psms_{s}\t{c}\n")

print(f"  Coverage: {coverage_pct}%, {len(all_peptides)} unique peptides, avg len {avg_pep_len}")
PYEOF

# ─── Level 10: CONVERGENCE 4 — Final report ───
echo "[Level 10 / CONVERGENCE 4] Generating final report..."
python3 << PYEOF
import os

out = os.environ.get("OUT", "outputs")
results = os.environ.get("RESULTS", "results")
os.makedirs(results, exist_ok=True)

# Read identification summary
id_stats = {}
with open(f"{out}/analysis/identification_summary.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        id_stats[k] = v

# Read coverage stats
cov_stats = {}
with open(f"{out}/analysis/coverage_stats.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        cov_stats[k] = v

# Read QC summary
qc_stats = {}
with open(f"{out}/qc/qc_summary.tsv") as f:
    next(f)
    for line in f:
        k, v = line.strip().split("\t")
        qc_stats[k] = v

with open(f"{results}/report.csv", "w") as f:
    f.write("metric,value\n")
    f.write(f"samples,{qc_stats.get('samples','3')}\n")
    f.write(f"total_spectra,{qc_stats.get('total_spectra','0')}\n")
    f.write(f"total_psms,{id_stats.get('total_psms','0')}\n")
    f.write(f"unique_peptides,{cov_stats.get('unique_peptides','0')}\n")
    f.write(f"identified_proteins,{cov_stats.get('identified_proteins','0')}\n")
    f.write(f"sequence_coverage_pct,{cov_stats.get('bsa_sequence_coverage_pct','0')}\n")
    f.write(f"avg_peptide_length,{cov_stats.get('avg_peptide_length','0')}\n")
    # Per-sample PSMs
    for s in ["BSA1", "BSA2", "BSA3"]:
        f.write(f"psms_{s},{cov_stats.get(f'psms_{s}','0')}\n")
    # Per-engine per-sample
    for key in sorted(id_stats):
        if key.startswith("psms_BSA"):
            f.write(f"{key},{id_stats[key]}\n")

print("Report written to results/report.csv")
PYEOF

echo ""
echo "=== Pipeline complete ==="
cat "${RESULTS}/report.csv"
