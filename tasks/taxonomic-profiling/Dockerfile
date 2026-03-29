FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

SHELL ["mamba", "run", "-n", "taxonomic-profiling", "/bin/bash", "-c"]

WORKDIR /work
COPY . /work/

RUN chmod +x run_script.sh
CMD ["mamba", "run", "-n", "taxonomic-profiling", "bash", "run_script.sh"]
