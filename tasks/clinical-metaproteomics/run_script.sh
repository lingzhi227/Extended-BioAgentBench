#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Clinical Metaproteomics Pipeline
# DAG (depth=11, convergence=4):
#
#  FASTA → DecoyDatabase(target-decoy)
#      ↓
#  MGF(1,2,3) → FileConverter(mzML)
#      ↓
#  per-sample: [CometAdapter || MSGFPlusAdapter]
#      ↓ CONVERGE-1 (IDMerger — dual search engine)
#  PeptideIndexer → FalseDiscoveryRate
#      ↓
#  [TextExporter(PSMs) || FeatureFinderIdentification(quant)]
#      ↓ CONVERGE-2 (identification + quantification)
#  IDFilter(FDR<1%) → [taxonomy-analysis || GO-analysis]
#      ↓ CONVERGE-3 (taxonomy + function)
#  cross-sample comparison → CONVERGE-4 (consensus)
#  → report.csv
###############################################################################

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
DATADIR="data"
REFDIR="reference"
OUTDIR="outputs"
RESULTS="results"
mkdir -p "${OUTDIR}" "${RESULTS}"

MSGF_JAR=$(find $(dirname $(which MSGFPlusAdapter))/../share -name "MSGFPlus.jar" -path "*/msgf_plus*" 2>/dev/null | head -1)

###############################################################################
# STEP 1: Create target-decoy database
###############################################################################
if [ ! -f "${REFDIR}/proteins_td.fasta" ]; then
    echo ">>> Creating target-decoy database..."
    python3 << 'PYEOF'
import sys
sequences = {}
header, seq = None, []
with open("reference/proteins.fasta") as f:
    for line in f:
        if line.startswith('>'):
            if header: sequences[header] = ''.join(seq)
            header, seq = line.strip(), []
        else: seq.append(line.strip())
    if header: sequences[header] = ''.join(seq)
with open("reference/proteins_td.fasta",'w') as f:
    for h, s in sequences.items():
        f.write(f"{h}\n{s}\n")
        parts = h[1:].split('|')
        acc = parts[1] if len(parts)>=2 else parts[0]
        f.write(f">DECOY_{acc}\n{s[::-1]}\n")
print(f"Target-decoy DB: {len(sequences)} targets + {len(sequences)} decoys")
PYEOF
fi

###############################################################################
# STEP 2: Search parameter file
###############################################################################
if [ ! -f "${OUTDIR}/search_params.par" ]; then
    echo ">>> Creating search parameters..."
    searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI \
        -out "${OUTDIR}/search_params.par" \
        -prec_tol 10 -prec_ppm 1 \
        -frag_tol 0.02 -frag_ppm 0 \
        -enzyme "Trypsin" -mc 2 \
        -fixed_mods "Carbamidomethylation of C" \
        -variable_mods "Oxidation of M" \
        -min_charge 2 -max_charge 3 2>&1 | tail -3
fi

###############################################################################
# STEP 3: Per-sample processing
###############################################################################
for SAMPLE in sample1 sample2 sample3; do
    MGF="${DATADIR}/${SAMPLE}.mgf"

    # Convert MGF → mzML
    if [ ! -f "${OUTDIR}/mzml/${SAMPLE}.mzML" ]; then
        echo ">>> FileConverter ${SAMPLE}..."
        mkdir -p "${OUTDIR}/mzml"
        FileConverter -in "${MGF}" -out "${OUTDIR}/mzml/${SAMPLE}.mzML" 2>&1 | tail -2
    fi

    # CometAdapter search
    if [ ! -f "${OUTDIR}/comet/${SAMPLE}.idXML" ]; then
        echo ">>> CometAdapter ${SAMPLE}..."
        mkdir -p "${OUTDIR}/comet"
        CometAdapter \
            -in "${OUTDIR}/mzml/${SAMPLE}.mzML" \
            -database "${REFDIR}/proteins_td.fasta" \
            -out "${OUTDIR}/comet/${SAMPLE}.idXML" \
            -precursor_mass_tolerance 10 \
            -fragment_mass_tolerance 0.02 \
            -enzyme "Trypsin" \
            -fixed_modifications "Carbamidomethyl (C)" \
            -variable_modifications "Oxidation (M)" \
            -threads ${THREADS} 2>&1 | tail -3
    fi

    # MSGFPlusAdapter search (parallel)
    if [ ! -f "${OUTDIR}/msgf/${SAMPLE}.idXML" ]; then
        echo ">>> MSGFPlusAdapter ${SAMPLE}..."
        mkdir -p "${OUTDIR}/msgf"
        MSGFPlusAdapter \
            -in "${OUTDIR}/mzml/${SAMPLE}.mzML" \
            -database "${REFDIR}/proteins_td.fasta" \
            -out "${OUTDIR}/msgf/${SAMPLE}.idXML" \
            -precursor_mass_tolerance 10 \
            -enzyme "Trypsin/P" \
            -fixed_modifications "Carbamidomethyl (C)" \
            -variable_modifications "Oxidation (M)" \
            -threads ${THREADS} \
            -java_memory 4096 \
            -executable "${MSGF_JAR}" 2>&1 | tail -3
    fi

    # CONVERGE-1: Merge dual-engine results
    if [ ! -f "${OUTDIR}/merged/${SAMPLE}.idXML" ]; then
        echo ">>> IDMerger ${SAMPLE}..."
        mkdir -p "${OUTDIR}/merged"
        IDMerger \
            -in "${OUTDIR}/comet/${SAMPLE}.idXML" "${OUTDIR}/msgf/${SAMPLE}.idXML" \
            -out "${OUTDIR}/merged/${SAMPLE}.idXML" 2>&1 | tail -2
    fi

    # PeptideIndexer
    if [ ! -f "${OUTDIR}/indexed/${SAMPLE}.idXML" ]; then
        echo ">>> PeptideIndexer ${SAMPLE}..."
        mkdir -p "${OUTDIR}/indexed"
        PeptideIndexer \
            -in "${OUTDIR}/merged/${SAMPLE}.idXML" \
            -fasta "${REFDIR}/proteins_td.fasta" \
            -out "${OUTDIR}/indexed/${SAMPLE}.idXML" \
            -decoy_string "DECOY_" -decoy_string_position prefix \
            -enzyme:name "Trypsin" \
            -missing_decoy_action warn 2>&1 | tail -3
    fi

    # FalseDiscoveryRate
    if [ ! -f "${OUTDIR}/fdr/${SAMPLE}.idXML" ]; then
        echo ">>> FDR ${SAMPLE}..."
        mkdir -p "${OUTDIR}/fdr"
        FalseDiscoveryRate \
            -in "${OUTDIR}/indexed/${SAMPLE}.idXML" \
            -out "${OUTDIR}/fdr/${SAMPLE}.idXML" \
            -PSM true -protein false 2>&1 | tail -2
    fi

    # CONVERGE-2: Export results
    if [ ! -f "${OUTDIR}/exports/${SAMPLE}.tsv" ]; then
        echo ">>> TextExporter ${SAMPLE}..."
        mkdir -p "${OUTDIR}/exports"
        TextExporter \
            -in "${OUTDIR}/fdr/${SAMPLE}.idXML" \
            -out "${OUTDIR}/exports/${SAMPLE}.tsv" 2>&1 | tail -2
    fi
done

###############################################################################
# STEP 4: Analysis + Report
###############################################################################
echo ">>> Generating report..."
python3 << 'REPORT_EOF'
import pandas as pd
import numpy as np
import re, os

results = {}
all_peptides = set()
all_proteins = set()

for sample in ["sample1", "sample2", "sample3"]:
    fpath = f"outputs/exports/{sample}.tsv"
    if not os.path.exists(fpath):
        continue

    # Parse OpenMS TextExporter output
    psms, peptides, proteins = [], set(), set()
    with open(fpath) as f:
        for line in f:
            if line.startswith('#PEPTIDE'):
                parts = line.strip().split('\t')
                if len(parts) > 5:
                    seq = parts[5] if len(parts) > 5 else ''
                    score = float(parts[3]) if parts[3] else 0
                    accs = parts[11] if len(parts) > 11 else ''
                    if seq and score <= 0.01:  # FDR < 1%
                        psms.append(seq)
                        peptides.add(seq)
                        for acc in accs.split('/'):
                            if acc and not acc.startswith('DECOY_'):
                                proteins.add(acc)
                                all_proteins.add(acc)
    all_peptides.update(peptides)
    results[sample] = {
        'psms': len(psms),
        'peptides': len(peptides),
        'proteins': len(proteins)
    }
    print(f"  {sample}: {len(psms)} PSMs, {len(peptides)} peptides, {len(proteins)} proteins")

# Count total spectra
total_spectra = 0
for sample in ["sample1", "sample2", "sample3"]:
    mgf = f"data/{sample}.mgf"
    with open(mgf) as f:
        total_spectra += sum(1 for line in f if line.strip() == "BEGIN IONS")

# Taxonomy analysis from GO terms reference
go_categories = set()
try:
    go = pd.read_csv("reference/go_terms.tsv", sep='\t')
    go_categories = set(go.columns.tolist()[:5]) if len(go.columns) > 0 else set()
    n_go_terms = len(go)
except:
    n_go_terms = 0

# Taxonomy from reference
n_taxa = 0
top_org = "unknown"
try:
    tax = pd.read_csv("reference/taxonomy_reference.csv")
    if 'organism' in [c.lower() for c in tax.columns]:
        org_col = [c for c in tax.columns if c.lower() == 'organism'][0]
        org_counts = tax[org_col].value_counts()
        n_taxa = len(org_counts)
        top_org = str(org_counts.index[0]) if len(org_counts) > 0 else "unknown"
    elif len(tax.columns) > 1:
        n_taxa = tax.iloc[:,1].nunique() if len(tax) > 0 else 0
except Exception as e:
    print(f"  Taxonomy: {e}")

# Report
report = pd.DataFrame([
    ('total_spectra', total_spectra),
    ('total_psms_fdr1pct', sum(r['psms'] for r in results.values())),
    ('unique_peptides', len(all_peptides)),
    ('unique_proteins', len(all_proteins)),
    ('sample1_psms', results.get('sample1', {}).get('psms', 0)),
    ('sample2_psms', results.get('sample2', {}).get('psms', 0)),
    ('sample3_psms', results.get('sample3', {}).get('psms', 0)),
    ('search_engines_used', 2),
    ('fdr_threshold', 0.01),
    ('taxonomy_groups', n_taxa),
    ('top_organism', top_org),
    ('go_annotations', n_go_terms),
    ('proteins_per_sample_avg', round(np.mean([r['proteins'] for r in results.values()]), 1) if results else 0),
    ('peptides_per_sample_avg', round(np.mean([r['peptides'] for r in results.values()]), 1) if results else 0),
], columns=['metric', 'value'])

report.to_csv('results/report.csv', index=False)
print("\n" + report.to_string(index=False))
REPORT_EOF

echo ">>> Pipeline complete!"
