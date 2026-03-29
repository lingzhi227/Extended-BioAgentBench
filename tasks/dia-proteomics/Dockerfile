FROM condaforge/mambaforge:latest
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy
SHELL ["mamba", "run", "-n", "dia-proteomics", "/bin/bash", "-c"]
WORKDIR /workspace
COPY . /workspace/
CMD ["mamba", "run", "-n", "dia-proteomics", "bash", "run_script.sh"]
