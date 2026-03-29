FROM condaforge/mambaforge:latest
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy
ENV PATH=/opt/conda/envs/rna-editing-detection/bin:$PATH
WORKDIR /workspace
COPY run_script.sh /workspace/
