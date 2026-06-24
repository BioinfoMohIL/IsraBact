# Lightweight image for parse_amr + plot_heatmap (pure Python).
# Phylogeny tools (Prokka, Roary, IQ-TREE) run in their own StaPH-B images,
# so this image only needs the Python stack.
# Build from the repo root (where this Dockerfile and scripts/ live):
#   docker build -t bioinfomoh/amr-heatmap:latest .
FROM python:3.11-slim-bookworm

RUN pip install --no-cache-dir pandas numpy matplotlib biopython

COPY scripts/ /opt/amr-heatmap/scripts/
RUN chmod +x /opt/amr-heatmap/scripts/*.py
ENV PATH=/opt/amr-heatmap/scripts:$PATH

RUN python -c "import pandas, numpy, matplotlib, Bio; print('py deps ok')"

WORKDIR /data
