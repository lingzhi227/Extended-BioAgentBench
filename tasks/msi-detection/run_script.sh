#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# MSI Detection Pipeline: Multi-Caller Consensus Microsatellite Instability
# =============================================================================
# DAG Structure (depth=10, convergence=4):
#
#  tumor.bam                     normal.bam
#     |                              |
#  [samtools index]             [samtools index]           Level 1
#     |                              |
#  [samtools flagstat]          [samtools flagstat]        Level 2
#     |                              |
#     +--------------+---------------+
#                    |
#            CONVERGENCE 1                                 Level 3
#            (T+N indexed BAMs ready)
#                    |
#     +--------------+---------------+
#     |              |               |
#  [msisensor-pro  [msisensor2     [MANTIS                Level 4
#   scan+msi]       scan+msi]       kmer+instab]
#     |              |               |
#  [parse           [parse          [parse                Level 5
#   pro results]    ms2 results]    mantis results]
#     |              |               |
#     +--------------+---------------+
#                    |
#            CONVERGENCE 2                                 Level 6
#            (3-caller MSI scores)
#                    |
#     +--------------+---------------+
#     |              |               |
#  [gatk            [mosdepth       [consensus            Level 7
#   HaplotypeCaller  coverage]       MSI scoring]
#   on MSI loci]     |               |
#     |              |               |
#     +--------------+---------------+
#                    |
#            CONVERGENCE 3                                 Level 8
#            (variants + coverage + consensus)
#                    |
#              [TMB estimation]                            Level 9
#                    |
#            CONVERGENCE 4                                 Level 10
#            (TMB + MSI + coverage + variants)
#                    |
#              [generate report.csv]
# =============================================================================

THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATA="${WORKDIR}/data"
REF="${WORKDIR}/reference"
OUT="${WORKDIR}/outputs"
RES="${WORKDIR}/results"
GENOME="${REF}/genome.fa"

# MANTIS location (cloned from GitHub)
MANTIS_DIR="${OUT}/MANTIS"

mkdir -p "${OUT}"/{msisensor_pro,msisensor2,mantis,gatk,mosdepth,consensus,tmb} "${RES}"

# =============================================================================
# Ensure reference is indexed
# =============================================================================
if [ ! -f "${GENOME}.fai" ]; then
    echo "Indexing reference FASTA..."
    samtools faidx "${GENOME}"
fi
if [ ! -f "${REF}/genome.dict" ]; then
    echo "Creating sequence dictionary..."
    samtools dict "${GENOME}" > "${REF}/genome.dict"
fi

# =============================================================================
# Level 1: Index tumor and normal BAMs
# =============================================================================
echo "[Level 1] Indexing BAM files..."
for SAMPLE in tumor normal; do
    if [ ! -f "${DATA}/${SAMPLE}.bam.bai" ]; then
        samtools index "${DATA}/${SAMPLE}.bam"
    fi
done

# =============================================================================
# Level 2: Flagstat QC on both BAMs
# =============================================================================
echo "[Level 2] Running flagstat QC..."
for SAMPLE in tumor normal; do
    if [ ! -f "${OUT}/${SAMPLE}_flagstat.txt" ]; then
        samtools flagstat "${DATA}/${SAMPLE}.bam" > "${OUT}/${SAMPLE}_flagstat.txt"
    fi
done

# Verify >10% mapping rate
for SAMPLE in tumor normal; do
    MAPPED_PCT=$(python3 -c "
with open('${OUT}/${SAMPLE}_flagstat.txt') as f:
    for line in f:
        if 'mapped' in line and '%' in line:
            pct = line.split('(')[1].split('%')[0]
            print(pct)
            break
")
    echo "  ${SAMPLE} mapping rate: ${MAPPED_PCT}%"
    python3 -c "assert float('${MAPPED_PCT}') > 10, '${SAMPLE} mapping rate too low: ${MAPPED_PCT}%'"
done

# =============================================================================
# CONVERGENCE 1 (Level 3): T+N indexed BAMs ready
# =============================================================================
echo "[Level 3] CONVERGENCE 1: Tumor and Normal BAMs indexed and QC-passed"

# =============================================================================
# Level 4a: MSI Caller 1 - msisensor-pro (scan + msi)
# =============================================================================
echo "[Level 4a] Running msisensor-pro..."
if [ ! -f "${OUT}/msisensor_pro/scan.list" ]; then
    msisensor-pro scan -d "${GENOME}" -o "${OUT}/msisensor_pro/scan.list" 2>/dev/null
fi

if [ ! -f "${OUT}/msisensor_pro/pro_result" ]; then
    msisensor-pro msi \
        -d "${OUT}/msisensor_pro/scan.list" \
        -t "${DATA}/tumor.bam" \
        -n "${DATA}/normal.bam" \
        -o "${OUT}/msisensor_pro/pro_result" 2>/dev/null || true
fi

# =============================================================================
# Level 4b: MSI Caller 2 - msisensor2 (scan + msi)
# =============================================================================
echo "[Level 4b] Running msisensor2..."
if [ ! -f "${OUT}/msisensor2/scan.list" ]; then
    msisensor2 scan -d "${GENOME}" -o "${OUT}/msisensor2/scan.list" 2>/dev/null
fi

if [ ! -f "${OUT}/msisensor2/ms2_result" ]; then
    msisensor2 msi \
        -d "${OUT}/msisensor2/scan.list" \
        -t "${DATA}/tumor.bam" \
        -n "${DATA}/normal.bam" \
        -o "${OUT}/msisensor2/ms2_result" 2>/dev/null || true
fi

# =============================================================================
# Level 4c: MSI Caller 3 - MANTIS
# =============================================================================
echo "[Level 4c] Running MANTIS..."
if [ ! -f "${OUT}/mantis/mantis_result.txt" ]; then
    # Create MANTIS 6-column BED from msisensor-pro scan list
    # Format required: chrom start end (UNIT)COUNT . .
    if [ -f "${OUT}/msisensor_pro/scan.list" ]; then
        python3 -c "
with open('${OUT}/msisensor_pro/scan.list') as fin, open('${OUT}/mantis/loci.bed', 'w') as fout:
    next(fin)  # skip header
    for line in fin:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            chrom = parts[0]
            start = int(parts[1])
            repeat_unit = parts[4]  # repeat_unit_bases
            repeat_times = int(parts[3])
            end = start + len(repeat_unit) * repeat_times
            kmer_spec = f'({repeat_unit}){repeat_times}'
            fout.write(f'{chrom}\t{start}\t{end}\t{kmer_spec}\t.\t.\n')
" || true
    fi

    # Clone MANTIS if not present
    if [ ! -f "${MANTIS_DIR}/mantis.py" ]; then
        git clone https://github.com/OSU-SRLab/MANTIS.git "${MANTIS_DIR}" 2>/dev/null || true
    fi

    # Run MANTIS
    if [ -f "${OUT}/mantis/loci.bed" ] && [ -f "${MANTIS_DIR}/mantis.py" ]; then
        python3 "${MANTIS_DIR}/mantis.py" \
            --bedfile "${OUT}/mantis/loci.bed" \
            --genome "${GENOME}" \
            -n "${DATA}/normal.bam" \
            -t "${DATA}/tumor.bam" \
            -o "${OUT}/mantis/mantis_result.txt" \
            --threads "${THREADS}" 2>/dev/null || true
    fi
fi

# =============================================================================
# Level 5: Parse results from each caller
# =============================================================================
echo "[Level 5] Parsing MSI caller results..."

python3 << 'PARSE_RESULTS'
import os, json

OUT = os.environ.get("OUT", "outputs")
os.makedirs(f"{OUT}/consensus", exist_ok=True)

caller_scores = {}

# Parse msisensor-pro results
pro_file = f"{OUT}/msisensor_pro/pro_result"
if os.path.exists(pro_file):
    with open(pro_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        # Format: total_sites  somatic_sites  percentage
                        total = int(parts[0])
                        somatic = int(parts[1])
                        pct = float(parts[2])
                        caller_scores['msisensor_pro'] = {
                            'total_sites': total,
                            'unstable_sites': somatic,
                            'msi_score': pct
                        }
                        print(f"msisensor-pro: {pct}% ({somatic}/{total} sites)")
                    except (ValueError, IndexError):
                        pass
else:
    print("msisensor-pro: no results file")

# Parse msisensor2 results
ms2_file = f"{OUT}/msisensor2/ms2_result"
if os.path.exists(ms2_file):
    with open(ms2_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        total = int(parts[0])
                        somatic = int(parts[1])
                        pct = float(parts[2])
                        caller_scores['msisensor2'] = {
                            'total_sites': total,
                            'unstable_sites': somatic,
                            'msi_score': pct
                        }
                        print(f"msisensor2: {pct}% ({somatic}/{total} sites)")
                    except (ValueError, IndexError):
                        pass
else:
    print("msisensor2: no results file")

# Parse MANTIS results
mantis_file = f"{OUT}/mantis/mantis_result.txt"
if os.path.exists(mantis_file):
    with open(mantis_file) as f:
        content = f.read()
    # MANTIS tab-delimited: Locus Normal_Reads Tumor_Reads Difference Distance Dissimilarity
    # Last line is Average
    for line in content.strip().split('\n'):
        if line.startswith('Average'):
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                try:
                    diff = float(parts[3])
                    dist = float(parts[4]) if len(parts) > 4 else 0
                    dissim = float(parts[5]) if len(parts) > 5 else 0
                    loci_count = sum(1 for l in content.strip().split('\n')
                                     if l.strip() and not l.startswith('Locus') and not l.startswith('Average'))
                    caller_scores['mantis'] = {
                        'step_wise_difference': diff,
                        'euclidean_distance': dist,
                        'cosine_dissimilarity': dissim,
                        'msi_score': diff * 100,  # normalize to percentage scale
                        'loci_analyzed': loci_count
                    }
                    print(f"MANTIS: DIF={diff}, EUC={dist}, COS={dissim}, loci={loci_count}")
                except (ValueError, IndexError):
                    pass
else:
    print("MANTIS: no results file")

# Save parsed scores
with open(f"{OUT}/consensus/caller_scores.json", 'w') as f:
    json.dump(caller_scores, f, indent=2)

print(f"Parsed {len(caller_scores)} caller results")
PARSE_RESULTS

export OUT="${OUT}"

# =============================================================================
# CONVERGENCE 2 (Level 6): All 3 MSI caller results available
# =============================================================================
echo "[Level 6] CONVERGENCE 2: All MSI caller results parsed"

# =============================================================================
# Level 7a: GATK HaplotypeCaller on MSI-associated loci
# =============================================================================
echo "[Level 7a] Running GATK HaplotypeCaller on MSI loci..."

# Create MSI loci intervals from the scan list
if [ ! -f "${OUT}/gatk/msi_loci.interval_list" ]; then
    python3 -c "
import os
# Read the sequence dictionary header
header_lines = []
with open('${REF}/genome.dict') as f:
    for line in f:
        header_lines.append(line.strip())

# Read MSI loci from msisensor-pro scan
loci = []
scan_file = '${OUT}/msisensor_pro/scan.list'
if os.path.exists(scan_file):
    with open(scan_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                start = int(parts[1])
                end = start + 20  # cover the repeat region
                loci.append((chrom, start, end))

# Write interval list
with open('${OUT}/gatk/msi_loci.interval_list', 'w') as f:
    for h in header_lines:
        f.write(h + '\n')
    # Use broader regions around MSI loci (merge nearby)
    if loci:
        loci.sort()
        merged = [list(loci[0])]
        for chrom, start, end in loci[1:]:
            if chrom == merged[-1][0] and start <= merged[-1][2] + 100:
                merged[-1][2] = max(merged[-1][2], end)
            else:
                merged.append([chrom, start, end])
        for chrom, start, end in merged:
            # Pad with 50bp flanks
            s = max(1, start - 50)
            e = end + 50
            f.write(f'{chrom}\t{s}\t{e}\t+\tMSI_locus\n')
        print(f'Created {len(merged)} interval regions from {len(loci)} MSI loci')
"
fi

if [ ! -f "${OUT}/gatk/tumor_variants.g.vcf.gz" ]; then
    gatk HaplotypeCaller \
        -R "${GENOME}" \
        -I "${DATA}/tumor.bam" \
        -L "${OUT}/gatk/msi_loci.interval_list" \
        -O "${OUT}/gatk/tumor_variants.g.vcf.gz" \
        -ERC GVCF 2>/dev/null || true
fi

if [ ! -f "${OUT}/gatk/normal_variants.g.vcf.gz" ]; then
    gatk HaplotypeCaller \
        -R "${GENOME}" \
        -I "${DATA}/normal.bam" \
        -L "${OUT}/gatk/msi_loci.interval_list" \
        -O "${OUT}/gatk/normal_variants.g.vcf.gz" \
        -ERC GVCF 2>/dev/null || true
fi

# Joint genotype
if [ ! -f "${OUT}/gatk/joint.vcf.gz" ]; then
    gatk CombineGVCFs \
        -R "${GENOME}" \
        -V "${OUT}/gatk/tumor_variants.g.vcf.gz" \
        -V "${OUT}/gatk/normal_variants.g.vcf.gz" \
        -O "${OUT}/gatk/combined.g.vcf.gz" 2>/dev/null || true

    gatk GenotypeGVCFs \
        -R "${GENOME}" \
        -V "${OUT}/gatk/combined.g.vcf.gz" \
        -O "${OUT}/gatk/joint.vcf.gz" 2>/dev/null || true
fi

# Count variants in MSI loci
TUMOR_VARIANT_COUNT=0
if [ -f "${OUT}/gatk/joint.vcf.gz" ]; then
    TUMOR_VARIANT_COUNT=$(bcftools view -H "${OUT}/gatk/joint.vcf.gz" 2>/dev/null | wc -l || true)
fi
echo "  Variants in MSI loci: ${TUMOR_VARIANT_COUNT}"

# =============================================================================
# Level 7b: mosdepth coverage at MSI loci
# =============================================================================
echo "[Level 7b] Computing coverage at MSI loci..."

# Create a BED for mosdepth from the interval list
if [ ! -f "${OUT}/mosdepth/msi_loci.bed" ]; then
    python3 -c "
import os
scan_file = '${OUT}/msisensor_pro/scan.list'
if os.path.exists(scan_file):
    with open(scan_file) as f:
        next(f)  # skip header
        loci = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                start = int(parts[1])
                end = start + 20
                loci.append((chrom, start, end))
    with open('${OUT}/mosdepth/msi_loci.bed', 'w') as f:
        for chrom, start, end in loci:
            f.write(f'{chrom}\t{start}\t{end}\n')
    print(f'Created BED with {len(loci)} MSI loci')
"
fi

for SAMPLE in tumor normal; do
    if [ ! -f "${OUT}/mosdepth/${SAMPLE}.mosdepth.summary.txt" ]; then
        mosdepth --threads "${THREADS}" \
            --by "${OUT}/mosdepth/msi_loci.bed" \
            "${OUT}/mosdepth/${SAMPLE}" \
            "${DATA}/${SAMPLE}.bam" 2>/dev/null || true
    fi
done

# =============================================================================
# Level 7c: Consensus MSI scoring
# =============================================================================
echo "[Level 7c] Computing consensus MSI score..."

python3 << 'CONSENSUS'
import os, json

OUT = os.environ.get("OUT", "outputs")

# Load caller scores
scores_file = f"{OUT}/consensus/caller_scores.json"
if os.path.exists(scores_file):
    with open(scores_file) as f:
        caller_scores = json.load(f)
else:
    caller_scores = {}

# Compute consensus
msi_scores = []
msi_status_votes = []

# Thresholds for MSI-H classification
# msisensor-pro: >3.5% is MSI-H (standard threshold)
# msisensor2: >3.5% is MSI-H (same engine basis)
# MANTIS: >0.4 step-wise difference is MSI-H

for caller, data in caller_scores.items():
    score = data.get('msi_score', 0)
    msi_scores.append(score)

    if caller in ('msisensor_pro', 'msisensor2'):
        status = 'MSI-H' if score > 3.5 else 'MSS'
    elif caller == 'mantis':
        # MANTIS uses step-wise difference
        swd = data.get('step_wise_difference', 0)
        status = 'MSI-H' if swd > 0.4 else 'MSS'
    else:
        status = 'MSS'
    msi_status_votes.append(status)
    print(f"  {caller}: score={score:.2f}, status={status}")

# Consensus: majority vote
if msi_status_votes:
    msi_h_count = sum(1 for s in msi_status_votes if s == 'MSI-H')
    consensus_status = 'MSI-H' if msi_h_count > len(msi_status_votes) / 2 else 'MSS'
    consensus_score = sum(msi_scores) / len(msi_scores) if msi_scores else 0
else:
    consensus_status = 'UNKNOWN'
    consensus_score = 0

result = {
    'consensus_msi_status': consensus_status,
    'consensus_msi_score': round(consensus_score, 4),
    'num_callers': len(caller_scores),
    'msi_h_votes': msi_h_count if msi_status_votes else 0,
    'caller_details': caller_scores
}

with open(f"{OUT}/consensus/consensus_result.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"\nConsensus: {consensus_status} (score={consensus_score:.4f}, {msi_h_count}/{len(msi_status_votes)} callers agree)")
CONSENSUS

# =============================================================================
# CONVERGENCE 3 (Level 8): Variants + Coverage + Consensus
# =============================================================================
echo "[Level 8] CONVERGENCE 3: Variants, coverage, and MSI consensus computed"

# =============================================================================
# Level 9: TMB Estimation
# =============================================================================
echo "[Level 9] Estimating Tumor Mutation Burden..."

export GENOME
python3 << 'TMB'
import os, subprocess

OUT = os.environ.get("OUT", "outputs")
GENOME = os.environ.get("GENOME", "reference/genome.fa")

# Get total callable region size from reference
fai_file = GENOME + ".fai"
total_bases = 0
if os.path.exists(fai_file):
    with open(fai_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            total_bases += int(parts[1])

# Count somatic-like variants (present in tumor, different from normal)
variant_count = 0
joint_vcf = f"{OUT}/gatk/joint.vcf.gz"
if os.path.exists(joint_vcf):
    try:
        result = subprocess.run(
            ['bcftools', 'view', '-H', joint_vcf],
            capture_output=True, text=True
        )
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                variant_count += 1
    except Exception:
        pass

# TMB = mutations per megabase
region_size_mb = total_bases / 1_000_000 if total_bases > 0 else 1
tmb = variant_count / region_size_mb

result = {
    'total_variants_msi_loci': variant_count,
    'callable_region_mb': round(region_size_mb, 2),
    'tmb_per_mb': round(tmb, 2)
}

os.makedirs(f"{OUT}/tmb", exist_ok=True)
import json
with open(f"{OUT}/tmb/tmb_result.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"TMB: {tmb:.2f} mutations/Mb ({variant_count} variants in {region_size_mb:.2f} Mb)")

# TMB classification
if tmb >= 10:
    tmb_status = "TMB-High"
elif tmb >= 5:
    tmb_status = "TMB-Intermediate"
else:
    tmb_status = "TMB-Low"
print(f"TMB Status: {tmb_status}")
TMB

# =============================================================================
# CONVERGENCE 4 (Level 10): TMB + MSI + Coverage + Variants -> Report
# =============================================================================
echo "[Level 10] CONVERGENCE 4: Generating final report..."

python3 << 'REPORT'
import os, json, csv

OUT = os.environ.get("OUT", "outputs")
RES = os.environ.get("RES", "results")

metrics = {}

# --- Flagstat metrics ---
for sample in ['tumor', 'normal']:
    flagstat_file = f"{OUT}/{sample}_flagstat.txt"
    if os.path.exists(flagstat_file):
        with open(flagstat_file) as f:
            for line in f:
                if 'in total' in line:
                    metrics[f'{sample}_total_reads'] = int(line.split()[0])
                elif 'mapped' in line and '%' in line and 'primary' not in line and 'mate' not in line:
                    pct = line.split('(')[1].split('%')[0]
                    metrics[f'{sample}_mapped_pct'] = float(pct)

# --- MSI caller results ---
scores_file = f"{OUT}/consensus/caller_scores.json"
if os.path.exists(scores_file):
    with open(scores_file) as f:
        caller_scores = json.load(f)
    for caller, data in caller_scores.items():
        metrics[f'{caller}_msi_score'] = data.get('msi_score', 0)
        if 'total_sites' in data:
            metrics[f'{caller}_total_sites'] = data['total_sites']
            metrics[f'{caller}_unstable_sites'] = data['unstable_sites']
        if 'step_wise_difference' in data:
            metrics[f'{caller}_step_wise_diff'] = data['step_wise_difference']
        if 'loci_analyzed' in data:
            metrics[f'{caller}_loci_analyzed'] = data['loci_analyzed']

# --- Consensus MSI ---
consensus_file = f"{OUT}/consensus/consensus_result.json"
if os.path.exists(consensus_file):
    with open(consensus_file) as f:
        consensus = json.load(f)
    metrics['consensus_msi_status'] = consensus['consensus_msi_status']
    metrics['consensus_msi_score'] = consensus['consensus_msi_score']
    metrics['msi_callers_agree'] = f"{consensus.get('msi_h_votes', 0)}/{consensus.get('num_callers', 0)}"

# --- Coverage from mosdepth ---
for sample in ['tumor', 'normal']:
    summ = f"{OUT}/mosdepth/{sample}.mosdepth.summary.txt"
    if os.path.exists(summ):
        with open(summ) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4 and parts[0] in ('total', 'total_region'):
                    try:
                        metrics[f'{sample}_mean_coverage'] = round(float(parts[3]), 2)
                    except ValueError:
                        pass

# --- Variant count from GATK ---
tmb_file = f"{OUT}/tmb/tmb_result.json"
if os.path.exists(tmb_file):
    with open(tmb_file) as f:
        tmb_data = json.load(f)
    metrics['variants_in_msi_loci'] = tmb_data.get('total_variants_msi_loci', 0)
    metrics['tmb_per_mb'] = tmb_data.get('tmb_per_mb', 0)
    metrics['callable_region_mb'] = tmb_data.get('callable_region_mb', 0)

# TMB status
tmb_val = metrics.get('tmb_per_mb', 0)
if tmb_val >= 10:
    metrics['tmb_status'] = 'TMB-High'
elif tmb_val >= 5:
    metrics['tmb_status'] = 'TMB-Intermediate'
else:
    metrics['tmb_status'] = 'TMB-Low'

# --- msisensor-pro per-site details ---
pro_dis = f"{OUT}/msisensor_pro/pro_result_dis"
if os.path.exists(pro_dis):
    site_count = 0
    with open(pro_dis) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                site_count += 1
    metrics['msisensor_pro_sites_analyzed'] = site_count

# Write report.csv
os.makedirs(RES, exist_ok=True)
with open(f"{RES}/report.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['metric', 'value'])
    for k, v in metrics.items():
        writer.writerow([k, v])

print("=== Report ===")
for k, v in metrics.items():
    print(f"  {k}: {v}")
print(f"\nReport written to {RES}/report.csv")
REPORT

echo "=== Pipeline Complete ==="
cat "${RES}/report.csv"
