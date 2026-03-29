FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml -n task && mamba clean -afy

SHELL ["mamba", "run", "-n", "task", "/bin/bash", "-c"]

COPY run_script.sh /app/run_script.sh
RUN chmod +x /app/run_script.sh

WORKDIR /app
CMD ["mamba", "run", "-n", "task", "bash", "run_script.sh"]
