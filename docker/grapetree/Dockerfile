FROM python:3.9-slim

LABEL maintainer="David Maimoun"
LABEL description="GrapeTree - Visualization of core genomic relationships"
LABEL version="2.2.0"

# 1. Installation des dépendances système
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# 2. Installation de pandas AVANT grapetree
# On ajoute pandas ici car grapetree en a absolument besoin pour s'importer
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    matplotlib \
    biopython \
    requests

# 3. Installation de GrapeTree
RUN pip install --no-cache-dir grapetree

# Create working directory
WORKDIR /data

# Test installation
RUN grapetree --help

# Default command
CMD ["grapetree", "--help"]