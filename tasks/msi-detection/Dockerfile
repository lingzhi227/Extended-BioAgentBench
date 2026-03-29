FROM condaforge/mambaforge:latest

WORKDIR /app

COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -afy

# Clone MANTIS for MSI detection
RUN git clone https://github.com/OSU-SRLab/MANTIS.git /app/MANTIS

SHELL ["mamba", "run", "-n", "msi-detection", "/bin/bash", "-c"]

COPY data/ data/
COPY reference/ reference/
COPY run_script.sh .

RUN chmod +x run_script.sh

CMD ["mamba", "run", "-n", "msi-detection", "bash", "run_script.sh"]
