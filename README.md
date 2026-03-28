# Extended-BioAgentBench

Extension of [BioAgentBench](https://github.com/bioagent-bench/bioagent-bench) with 10 additional bioinformatics benchmark tasks covering new domains and increasing workflow complexity.

# Contents

- **src/**
  - **dataset.py**: Code for downloading the input and truth files
- **tasks/**: Individual bioinformatics tasks
  - **task_name/**
  - **Dockerfile**: You can reproduce the truth files with this
  - **run_script.sh**: Script to run the bioinformatics tools
  - **environment.yml**: Conda environment specification

# Tasks

| Task ID | Domain |
|---------|--------|
| chipseq-peak-calling | Epigenomics |
| bacterial-assembly | Microbiology |
| mobile-elements | Microbiology |
| outbreak-investigation | Epidemiology |
| atacseq-accessibility | Epigenomics |
| longread-assembly | Microbiology |
| hybrid-assembly | Microbiology |
| sv-detection | Genomics |
| pangenome-evolution | Microbiology |
| metagenomic-profiling | Metagenomics |

# General instructions
### CLI: dataset downloader
Use the Click-based CLI in `src/dataset.py` to list tasks and download data, reference files, and results.

- **Install uv (if not installed)**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

- **Create a virtual environment and install dependencies**
```bash
uv sync
```

- **Show help**
```bash
uv run python src/dataset.py --help
```

- **Download all input data to a destination**
```bash
uv run python src/dataset.py download --all --dest /path/to/output/
```

- **Download all input data, reference data, results data to a destination**
```bash
uv run python src/dataset.py download --all --dest /path/to/output/ --reference --results
```

Other CLI options (see `--help` for details): list available tasks, download a single task, include reference files, limit output to results only, and combine multiple tasks in one call.

Files are downloaded under `tasks/<task_id>/`:
```
tasks/<task_id>/
  data/         # Input data (FASTQ, BED, etc.)
  reference/    # Reference genomes (if needed)
  results/      # Ground truth output
```

# Citation

Based on BioAgentBench: https://arxiv.org/abs/2601.21800
