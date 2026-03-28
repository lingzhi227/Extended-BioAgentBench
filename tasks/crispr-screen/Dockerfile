FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH /opt/conda/envs/crispr-screen/bin:$PATH
SHELL ["mamba", "run", "-n", "crispr-screen", "/bin/bash", "-c"]

WORKDIR /workspace
COPY run_script.sh /workspace/
RUN chmod +x /workspace/run_script.sh
