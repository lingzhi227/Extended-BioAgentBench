FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH /opt/conda/envs/dda-lfq-proteomics/bin:$PATH
SHELL ["conda", "run", "-n", "dda-lfq-proteomics", "/bin/bash", "-c"]

WORKDIR /workspace
COPY . /workspace/

CMD ["bash", "run_script.sh"]
