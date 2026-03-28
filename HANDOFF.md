# Extended-BioAgentBench Task Generation Handoff

## Context
Building a bioinformatics benchmark with 200 tasks total. Currently at 24 tasks (11-34).
Tasks 21, 31, 33, 34 are being fixed (replacing simulated data with real data).

## Your assignment
Generate benchmark tasks following the skill at:
`/global/homes/l/lingzhi/.claude/skills/bench-task-discovery/SKILL.md`

## Key paths
- Project: `/pscratch/sd/l/lingzhi/Extended-BioAgentBench/`
- Build dir: `/pscratch/sd/l/lingzhi/bench-task-output/`
- Metadata: `/pscratch/sd/l/lingzhi/Extended-BioAgentBench/src/task_metadata.json`
- Shared env: `/pscratch/sd/l/lingzhi/micromamba/envs/bench-all` (has most tools)
- HF token: `READ_FROM_KEY_FILE`
- OpenAI key: in `/pscratch/sd/l/lingzhi/key.md`
- GitHub: `lingzhi227/Extended-BioAgentBench`
- HuggingFace: `lingzhi227/Extended-BioAgentBench`

## Per-task lifecycle (MANDATORY, every step)
1. Create pluggable env: `micromamba create -p /pscratch/sd/l/lingzhi/micromamba/envs/task{N} ...`
2. Verify ALL tools: `which tool && tool --version`
3. Download REAL data (NEVER simulated/wgsim). Verify: size > 1KB, `head -1` is FASTQ/FASTA
4. Write run_script.sh with DAG comments and skip-if-exists guards
5. Run pipeline to completion. Debug ALL errors
6. Inspect results: all values present, no empty/N/A, no formatting issues
7. Package: `tar czf data.tar.gz -C data .` etc
8. Register in metadata (example FROM actual ground truth, leak check)
9. Upload to HuggingFace (archives + scripts + metadata)
10. Update README (task count, table, badge) + push to GitHub
11. Upload README to HuggingFace
12. Delete pluggable env

## Critical rules
- NEVER use simulated data (wgsim, art_illumina, etc)
- NEVER leak tool names in prompts
- Example values MUST come from actual pipeline output
- Comment DAG structure in run_script.sh
- Use `|| true` for non-critical grep/awk, NOT `|| echo "0"` (causes stray lines)
- MEGAHIT refuses to run if output dir exists (even empty) — delete before re-run
- Picard needs @RG header in BAM
- `grep -c` returns exit 1 when count=0 (breaks set -e). Use `grep -c ... || true; VAR=${VAR:-0}`

## Domains already covered (DO NOT duplicate)
ChIP-seq, ATAC-seq, bacterial assembly, mobile elements, outbreak investigation,
long-read assembly, hybrid assembly, SV detection, pan-genome, metagenomics,
phage, genome comparison, mapping QC, multi-sample variants, consensus genome,
gene prediction, downsampling, plasmid typing, genome completeness, species ID,
viral amplicon, bisulfite methylation, RNA-seq isoform, ancient DNA
