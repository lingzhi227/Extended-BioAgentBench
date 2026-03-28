# Extended-BioAgentBench Task Generation — Session Assignment

You are generating bioinformatics benchmark tasks for Extended-BioAgentBench.

## READ THESE FILES FIRST
1. Skill: `/global/homes/l/lingzhi/.claude/skills/bench-task-discovery/SKILL.md`
2. Existing metadata: `/pscratch/sd/l/lingzhi/Extended-BioAgentBench/src/task_metadata.json`
3. This prompt (contains critical mistake catalog)

## YOUR ASSIGNMENT
**Session {A/B/C/D/E}**: Generate tasks {N1} and {N2} in domains {domain1} and {domain2}.

Work in YOUR dedicated build directory:
- Session A: `/pscratch/sd/l/lingzhi/bench-task-output/session-a/`
- Session B: `/pscratch/sd/l/lingzhi/bench-task-output/session-b/`
- Session C: `/pscratch/sd/l/lingzhi/bench-task-output/session-c/`
- Session D: `/pscratch/sd/l/lingzhi/bench-task-output/session-d/`
- Session E: `/pscratch/sd/l/lingzhi/bench-task-output/session-e/`

Use YOUR dedicated pluggable envs:
- Session A: `/pscratch/sd/l/lingzhi/micromamba/envs/sessA-task{N}`
- Session B: `/pscratch/sd/l/lingzhi/micromamba/envs/sessB-task{N}`
- etc.

**DO NOT** touch other sessions' directories or envs.

## KEY PATHS (shared, read-only except metadata)
- Extended-BioAgentBench repo: `/pscratch/sd/l/lingzhi/Extended-BioAgentBench/`
- Metadata (append-only): `/pscratch/sd/l/lingzhi/Extended-BioAgentBench/src/task_metadata.json`
- Shared env (read-only, has many tools): `/pscratch/sd/l/lingzhi/micromamba/envs/bench-all`
- HF token: `READ_FROM_KEY_FILE`
- HF repo: `lingzhi227/Extended-BioAgentBench`
- GitHub repo: `lingzhi227/Extended-BioAgentBench`

## MANDATORY PER-TASK LIFECYCLE (DO EVERY STEP, DO NOT SKIP)

### Step 1: Create pluggable env
```bash
micromamba create -p /pscratch/sd/l/lingzhi/micromamba/envs/sess{X}-task{N} \
  -c bioconda -c conda-forge <tools> -y
```

### Step 2: Verify EVERY tool
```bash
for tool in <list>; do
  WHICH=$(micromamba run -p .../sess{X}-task{N} which $tool 2>/dev/null)
  echo "${tool}: ${WHICH:-NOT FOUND}"
done
```
**If any tool is NOT FOUND, DO NOT PROCEED. Fix first.**

### Step 3: Download REAL data + verify
```bash
curl -sL "<url>" -o data/reads_R1.fastq.gz
# MANDATORY VERIFICATION:
ls -lh data/reads_R1.fastq.gz  # MUST be > 1KB
zcat data/reads_R1.fastq.gz | head -1  # MUST start with @
file data/reads_R1.fastq.gz  # MUST say "gzip compressed", NOT "HTML"
```

### Step 4: Write run_script.sh
- Include DAG structure as comments at the top
- Use `set -euo pipefail`
- Use skip-if-exists guards
- Use `THREADS=$(( $(nproc) > 8 ? 8 : $(nproc) ))`

### Step 5: Run pipeline + debug
Run to completion. If it fails, FIX and re-run. Do NOT skip broken steps.

### Step 6: Inspect results
```bash
while IFS=, read -r m v; do
  [ "$m" != "metric" ] && echo "${m}=${v} [$([ -z "$v" ] && echo EMPTY || echo OK)]"
done < results/report.csv
grep -cE ",,|,$|N/A" results/report.csv  # MUST be 0
```

### Step 7: Package
```bash
tar czf data.tar.gz -h -C data .  # -h dereferences symlinks!
tar czf reference.tar.gz -h -C reference .
tar czf results.tar.gz -C results .
```

### Step 8: Register in metadata
```python
# Example MUST come from actual ground truth, not invented:
gt = open("results/report.csv").read()
new_task["task_prompt"] = "..." + "\n<example>" + gt + "</example>"
# Run leak check IMMEDIATELY after:
leak_words = [<full list>]
leaks = [w for w in leak_words if re.search(r'\b'+re.escape(w)+r'\b', prompt, re.I)]
assert not leaks
```

### Step 9: Upload to HuggingFace
```python
for f in ["data.tar.gz", "reference.tar.gz", "results.tar.gz"]:
    api.upload_file(path_or_fileobj=f, path_in_repo=f"tasks/{tid}/{f}", ...)
for f in ["run_script.sh", "environment.yml", "Dockerfile"]:
    api.upload_file(...)
api.upload_file(path_or_fileobj="src/task_metadata.json", ...)
```

### Step 10: Update README + push GitHub + push HF README
```python
# Update task table in README.md
# Update badge count
# git add -A && git commit && git push
# Upload README.md to HuggingFace
```

### Step 11: Delete pluggable env
```bash
micromamba env remove -p .../sess{X}-task{N} -y
```

---

## MISTAKE CATALOG (REAL EXAMPLES — DO NOT REPEAT)

### M1: SIMULATED DATA
**What happened**: Used `wgsim` to generate fake FASTQ reads from reference genomes.
**Why it's wrong**: Benchmark is for real-world use. Simulated reads have uniform quality scores (all `8` or `I`), no adapter contamination, no sequencing errors — unrealistic.
**How to avoid**: ALWAYS download from SRA/ENA/Zenodo. Verify with `zcat | head -1` — real reads have variable quality, real headers like `@SRR...` or `@ERR...`.

### M2: PROMPT LEAKS TOOL NAMES
**What happened**: Task prompt said "gene annotation with Prokka, IS element detection with ISEScan, AMR detection with AMRFinderPlus"
**Why it's wrong**: Tells the agent exactly which tools to use — no challenge.
**How to avoid**: Prompt says WHAT to produce, never HOW. Compare:
- BAD: "Use Prokka for annotation, BUSCO for completeness"
- GOOD: "Annotate the genome and assess completeness"

### M3: EXAMPLE VALUES DON'T MATCH GROUND TRUTH
**What happened**: Wrote fake example values in prompts BEFORE running the pipeline, then forgot to update them.
**Result**: Prompt said `n50=4700000` but ground truth was `n50=123189`.
**How to avoid**: ALWAYS write examples AFTER pipeline runs, by reading actual results file.

### M4: METRIC NAMES LEAK TOOLS
**What happened**: Ground truth CSV had `busco_completeness`, `mlst_sequence_type`, `checkv_quality`.
**Why it's wrong**: Metric names reveal which tools to use.
**How to avoid**: Use generic names: `completeness`, `sequence_type`, `quality_tier`.

### M5: `grep -c` EXIT CODE WITH set -e
**What happened**: `grep -c "pattern" file || echo "0"` prints TWO lines when count=0.
`grep -c` returns exit code 1 when count is 0, which triggers the `||` branch, printing "0" to stdout AFTER grep already printed "0".
**How to avoid**: Use `VAR=$(grep -c ... 2>/dev/null || true); VAR=${VAR:-0}`

### M6: MEGAHIT OUTPUT DIR
**What happened**: MEGAHIT refuses to run if output directory already exists, even empty.
**How to avoid**: Always `rm -rf outputs/assembly` before running MEGAHIT. Or use `--force` flag if available.

### M7: PICARD NEEDS @RG HEADER
**What happened**: `picard MarkDuplicates` crashes with NullPointerException if BAM has no read group.
**How to avoid**: Always add `--rg-id sample --rg "SM:sample"` to bowtie2, or `-R "@RG\tID:sample\tSM:sample"` to bwa mem.

### M8: ZENODO/ENA DOWNLOADS RETURN HTML
**What happened**: `curl -sL <url> -o file.fastq.gz` downloaded an HTML error page (14 bytes or 14KB).
**How to avoid**: After every download, check: `file data.fastq.gz` must say "gzip", NOT "HTML". Size must be > 1KB.

### M9: WRONG REFERENCE FOR DATA
**What happened**: Cave bear WGS reads aligned to mitogenome → 0% mapping. Phage reads from mixed sample → 504 fragmented contigs.
**How to avoid**: Verify reference matches data. After alignment, check `samtools flagstat` — if <10% mapping, something is wrong.

### M10: FTP URLS BLOCKED
**What happened**: `curl -sL "ftp://ftp.sra.ebi.ac.uk/..."` returns exit 78 on HPC.
**How to avoid**: Use HTTPS: `https://ftp.sra.ebi.ac.uk/...`

### M11: GLOB IN DOUBLE QUOTES
**What happened**: `REPORT="${OUT}/aligned/"*PE_report.txt` doesn't expand the glob.
**How to avoid**: Use `REPORT=$(ls ${OUT}/aligned/*PE_report.txt 2>/dev/null | head -1)`

### M12: TAB-DELIMITED PARSING
**What happened**: `awk -F: '{print $2}'` failed on tab-delimited output (bismark, samtools stats).
**How to avoid**: Check actual output format first. Use `cut -f3` or `awk -F'\t' '{print $NF}'`.

### M13: FORGOT TO UPDATE README + HUGGINGFACE
**What happened**: Added 12 new tasks but README still showed 10. HuggingFace had old README.
**How to avoid**: ALWAYS update README task table + badge count + push to BOTH GitHub and HuggingFace.

---

## DOMAINS ALREADY COVERED (DO NOT DUPLICATE)
ChIP-seq, ATAC-seq, bacterial assembly, mobile elements, outbreak investigation,
long-read assembly, hybrid assembly, SV detection, pan-genome, metagenomics,
phage characterization, genome comparison, mapping QC, multi-sample variants,
consensus genome, gene prediction, downsampling analysis, plasmid typing,
genome completeness, species identification, viral amplicon, bisulfite methylation,
RNA-seq isoform, ancient DNA

## WHEN YOU'RE DONE
Run the full verification:
```python
import json
with open("/pscratch/sd/l/lingzhi/Extended-BioAgentBench/src/task_metadata.json") as f:
    tasks = json.load(f)
for t in tasks:
    # Check prompt leaks, URLs, example match
    ...
```
Then confirm: README updated on GitHub AND HuggingFace, all tar.gz uploaded, env deleted.
