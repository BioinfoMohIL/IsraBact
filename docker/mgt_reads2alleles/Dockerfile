FROM continuumio/miniconda3:latest

WORKDIR /app

RUN apt-get update && apt-get install -y \
    git \
    wget \
    unzip \
    tar \
    && apt-get clean

RUN mkdir -p /app/db && \
    wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz -O /app/db/minikraken.tgz && \
    tar -xzf /app/db/minikraken.tgz -C /app/db && \
    rm /app/db/minikraken.tgz

RUN git clone https://github.com/LanLab/MGT_reads2alleles.git

ENV KRAKEN_DEFAULT_DB=/app/db

RUN conda update -n base -c defaults conda && \
    conda env create -f MGT_reads2alleles/fq_to_allele.yml -n fq_to_allele

RUN chmod 755 MGT_reads2alleles/shovill_cmd/bin/shovill_15cov

RUN echo '#!/bin/bash\nconda run -n fq_to_allele python /app/MGT_reads2alleles/reads_to_alleles.py "$@"' > /bin/reads_to_alleles && \
    chmod +x /bin/reads_to_alleles

RUN apt-get update && apt-get install -y file

ENTRYPOINT ["/bin/bash"]
