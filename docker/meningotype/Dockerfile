FROM condaforge/miniforge3:latest

LABEL maintainer="bioinfomoh"
LABEL description="Neisseria meningitidis typing (meningotype, py3 patched)"
LABEL version="1.0.3"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*

# Source officielle + environnement testé par leur CI
RUN git clone --depth 1 https://github.com/MDU-PHL/meningotype.git /opt/meningotype_src

RUN mamba env create -n meningotype -f /opt/meningotype_src/environment.yml \
    && mamba clean -afy

ENV PATH="/opt/conda/envs/meningotype/bin:${PATH}"

# Correctifs Python 3
COPY patch_meningotype.py /tmp/patch_meningotype.py
RUN python /tmp/patch_meningotype.py /opt/meningotype_src/meningotype

RUN pip install --no-cache-dir /opt/meningotype_src

RUN meningotype --checkdeps
RUN meningotype --updatedb
RUN meningotype --version \
    && meningotype --test --all

WORKDIR /data
CMD ["meningotype", "--help"]