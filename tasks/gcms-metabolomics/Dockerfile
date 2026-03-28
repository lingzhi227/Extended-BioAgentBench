FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH /opt/conda/envs/gcms-metabolomics/bin:$PATH
SHELL ["conda", "run", "-n", "gcms-metabolomics", "/bin/bash", "-c"]

WORKDIR /workspace
COPY . /workspace/

CMD ["bash", "run_script.sh"]
