FROM condaforge/mambaforge:latest
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy
SHELL ["mamba", "run", "-n", "longread-rna-isoform", "/bin/bash", "-c"]
WORKDIR /workspace
COPY . /workspace/
CMD ["mamba", "run", "-n", "longread-rna-isoform", "bash", "run_script.sh"]
