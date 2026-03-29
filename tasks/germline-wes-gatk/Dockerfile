FROM condaforge/mambaforge:latest
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy
ENV PATH=/opt/conda/envs/germline-wes-gatk/bin:$PATH
WORKDIR /workspace
COPY run_script.sh /workspace/
