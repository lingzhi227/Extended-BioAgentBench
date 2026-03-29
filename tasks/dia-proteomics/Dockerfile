FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH /opt/conda/envs/dia-proteomics/bin:$PATH
SHELL ["conda", "run", "-n", "dia-proteomics", "/bin/bash", "-c"]

WORKDIR /workspace
COPY . /workspace/

CMD ["bash", "run_script.sh"]
