FROM ubuntu:25.04

WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
    python3 python3-pip python3-dev \
    make gcc g++ \
    libxml2-dev libxslt1-dev \
    ncbi-blast+ \
    prodigal \
    git wget curl unzip && \
    apt-get clean

# Install chewBBACA
RUN pip install chewBBACA==2.8.5 \
    numpy \
    pandas \
    matplotlib \
    networkx \
    scipy \
    --break-system-packages

