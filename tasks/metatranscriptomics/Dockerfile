FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

SHELL ["mamba", "run", "-n", "metatranscriptomics", "/bin/bash", "-c"]

WORKDIR /workspace
COPY . /workspace/

CMD ["mamba", "run", "-n", "metatranscriptomics", "bash", "run_script.sh"]
