FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH /opt/conda/envs/mag-recovery/bin:$PATH
SHELL ["mamba", "run", "-n", "mag-recovery", "/bin/bash", "-c"]

WORKDIR /workspace
COPY run_script.sh /workspace/
RUN chmod +x /workspace/run_script.sh
