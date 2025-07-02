FROM ubuntu:22.04

LABEL maintainer="bioinfomoh@gmail.com"

# Install core system dependencies
RUN apt-get update && apt-get install -y \
    build-essential git wget curl unzip bc gawk \
    libssl-dev default-jdk python3 python3-pip python3-dev \
    autoconf automake libboost-all-dev zlib1g-dev g++ maven flash \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
ENV PATH="/opt/conda/bin:$PATH"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    conda init

# Create and activate conda environment, install tools in one step
RUN conda create -n bioinfo_env -c conda-forge -c bioconda -c defaults \
    python=3.10 prokka platon spades minimap2 bwa-mem2 \
    samtools bedtools diamond blast hmmer infernal gnuplot aragorn && \
    conda clean -a -y

ENV PATH="/opt/conda/envs/bioinfo_env/bin:$PATH"

# Install Pasteur GitLab tools

WORKDIR /opt

# Installation de FASTA2AGP
RUN git clone https://gitlab.pasteur.fr/GIPhy/FASTA2AGP.git && \
    cd FASTA2AGP/src && \
    javac FASTA2AGP.java && \
    echo "Main-Class: FASTA2AGP" > MANIFEST.MF && \
    jar -cmvf MANIFEST.MF FASTA2AGP.jar FASTA2AGP.class && \
    rm MANIFEST.MF FASTA2AGP.class && \
    printf '#!/bin/bash\nexec java -jar /opt/FASTA2AGP/src/FASTA2AGP.jar "$@"' > /bin/FASTA2AGP && \
    chmod +x /bin/FASTA2AGP


# Install SAM2MSA tools
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/GIPhy/SAM2MSA.git && \
    cd SAM2MSA/src && \
    for p in SAM2MAP MAP2FASTA FASTA2MSA; do \
      javac "$p.java" && \
      echo "Main-Class: $p" > MANIFEST.MF && \
      jar -cmvf MANIFEST.MF "$p.jar" "$p.class" && \
      rm MANIFEST.MF "$p.class"; \
      printf '#!/bin/bash\nexec java -jar /opt/SAM2MSA/src/%s.jar "$@"\n' "$p" > "/bin/$p" && \
      chmod +x "/bin/$p"; \
    done


## contig_info
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/GIPhy/contig_info.git && \
    install -m 755 contig_info/contig_info.sh /bin/contig_info

# fqCleanER
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/GIPhy/fqCleanER.git && \
    install -m 755 fqCleanER/fqCleanER.sh /bin/fqCleanER


# Install fqtools (fqstats)
WORKDIR /opt
RUN wget ftp://ftp.pasteur.fr/pub/gensoft/projects/fqtools/fqtools-1.1a.tar.gz && \
    tar -xzf fqtools-1.1a.tar.gz && \
    cd fqtools-1.1a && ./configure && make && make install && \
    cd .. && rm -rf fqtools-1.1a*

# Install DSRC
WORKDIR /opt
RUN git clone https://github.com/refresh-bio/DSRC.git && \
    cd DSRC && sed -i '/^rm dsrc$/d' src/Makefile && \
    make -f Makefile bin && \
    cp bin/dsrc /usr/local/bin/ && \
    cd .. && rm -rf DSRC

# Alien tools directory
RUN mkdir -p /bin/alien

# AlienDiscover
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/GIPhy/AlienDiscover.git && \
    cd AlienDiscover/src && \
    javac AlienDiscover.java && \
    echo "Main-Class: AlienDiscover" > MANIFEST.MF && \
    jar -cmvf MANIFEST.MF AlienDiscover.jar AlienDiscover.class && \
    mkdir -p /bin/alien && \
    mv AlienDiscover.jar /bin/alien/ && \
    printf '#!/bin/bash\nexec java -jar /bin/alien/AlienDiscover.jar "$@"\n' > /bin/AlienDiscover && \
    chmod +x /bin/AlienDiscover && \
    cd /opt && rm -rf AlienDiscover

# AlienRemover
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/GIPhy/AlienRemover.git && \
    cd AlienRemover/src && \
    javac LHBF.java AlienRemover.java && \
    echo "Main-Class: AlienRemover" > MANIFEST.MF && \
    jar -cmvf MANIFEST.MF AlienRemover.jar AlienRemover.class LHBF.class && \
    mv AlienRemover.jar /bin/alien/ && \
    printf '#!/bin/bash\nexec java -jar /bin/alien/AlienRemover.jar "$@"\n' > /bin/AlienRemover && \
    chmod +x /bin/AlienRemover && \
    cd /opt && rm -rf AlienRemover

# AlienTrimmer
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/GIPhy/AlienTrimmer.git && \
    cd AlienTrimmer/src && \
    javac AlienTrimmer.java && \
    echo "Main-Class: AlienTrimmer" > MANIFEST.MF && \
    jar -cmvf MANIFEST.MF AlienTrimmer.jar AlienTrimmer.class && \
    mv AlienTrimmer.jar /bin/alien/ && \
    printf '#!/bin/bash\nexec java -jar /bin/alien/AlienTrimmer.jar "$@"\n' > /bin/AlienTrimmer && \
    chmod +x /bin/AlienTrimmer && \
    cd /opt && rm -rf AlienTrimmer


# Musket
COPY tools/musket-1.1.tar.gz /opt/
WORKDIR /opt
RUN tar -xzf musket-1.1.tar.gz && \
    cd musket-1.1 && make && cp musket /bin/ && cd .. && rm -rf musket-1.1 musket-1.1.tar.gz

# ROCK
WORKDIR /opt
RUN git clone https://gitlab.pasteur.fr/vlegrand/ROCK.git && \
    cd ROCK && \
    autoreconf -fi && \
    ./configure --prefix=/usr/local CFLAGS='-O3' CXXFLAGS='-O3' && \
    make -j$(nproc) && \
    make check || echo "Warning: some tests failed, continuing anyway" && \
    make install && \
    cd .. && rm -rf ROCK

# Clone and build ntCard
WORKDIR /opt
RUN git clone https://github.com/bcgsc/ntCard.git && \
    cd ntCard && \
    ./autogen.sh && \
    ./configure CFLAGS='-g -O3' CXXFLAGS='-g -O3' && \
    make && \
    cp ntcard /bin/ && \
    cd .. && rm -rf ntCard

# Final setup
WORKDIR /app
COPY COPYING README.md fq2dna.sh /app/
RUN chmod +x /app/fq2dna.sh

COPY ./bin/fq2dna /bin
RUN chmod +x /bin/fq2dna


CMD ["bash"]
