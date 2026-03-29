FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH="/opt/conda/envs/immune-repertoire/bin:$PATH"

WORKDIR /workspace
COPY run_script.sh /workspace/
RUN chmod +x /workspace/run_script.sh
