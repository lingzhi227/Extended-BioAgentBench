FROM mambaorg/micromamba:1.5.1
COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml -y && micromamba clean -a -y
ENV PATH="/opt/conda/envs/env/bin:$PATH"
WORKDIR /workspace
COPY run_script.sh /workspace/
