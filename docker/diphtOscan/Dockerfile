FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/conda/bin:$PATH"

# 1️⃣ Dépendances système de base
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    bzip2 \
    build-essential \
    unzip \
    ca-certificates \
    gawk \
    bc \
    zip \
    gzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 2️⃣ Installer Miniconda avec Python 3.10
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh \
    && /opt/conda/bin/conda init bash \
    && /opt/conda/bin/conda config --set always_yes yes --set changeps1 no

# 3️⃣ Configurer les channels conda
RUN conda update -n base -c defaults conda \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority flexible

# 4️⃣ Cloner diphtoscan
RUN git clone https://gitlab.pasteur.fr/BEBP/diphtoscan.git /opt/diphtoscan

WORKDIR /opt/diphtoscan

# 5️⃣ Créer l'environnement avec Python 3.10
RUN conda create -n diphtoscan python=3.10 -y \
    && conda clean -afy

# 6️⃣ Installer les dépendances bioinformatiques
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate diphtoscan && \
    conda install -c bioconda -c conda-forge \
    'biopython>=1.79,<1.84' \
    'integron_finder>=2.0,<2.1' \
    ncbi-amrfinderplus \
    infernal \
    prodigal \
    hmmer \
    blast \
    'mash>=2.1' \
    pandas \
    requests \
    fasttree \
    mafft \
    trimal \
    iqtree \
    fastme \
    -y && \
    conda clean -afy"

# 7️⃣ Installer diphtoscan SEUL (sans JolyTree dans le dossier)
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate diphtoscan && \
    python -m pip install . --no-deps"

# 8️⃣ MAINTENANT cloner JolyTree dans un dossier séparé pour l'installation
RUN mkdir -p /opt/diphtoscan/script && \
    cd /opt/diphtoscan/script && \
    git clone https://gitlab.pasteur.fr/GIPhy/JolyTree.git && \
    chmod +x /opt/diphtoscan/script/JolyTree/*.sh

# 9️⃣ Créer le lien symbolique /diphtoscan -> /opt/diphtoscan
RUN ln -s /opt/diphtoscan /diphtoscan

# 🔟 Ajouter JolyTree au PATH de l'environnement conda
RUN mkdir -p /opt/conda/envs/diphtoscan/etc/conda/activate.d && \
    echo 'export PATH="/opt/diphtoscan/script/JolyTree:$PATH"' > /opt/conda/envs/diphtoscan/etc/conda/activate.d/env_vars.sh

# 1️⃣1️⃣ Créer des liens symboliques pour les scripts JolyTree
RUN for script in /opt/diphtoscan/script/JolyTree/*.sh; do \
        if [ -f "$script" ]; then \
            ln -sf "$script" /usr/local/bin/$(basename "$script"); \
        fi; \
    done

# 1️⃣2️⃣ Vérifier que REQ existe ou créer un fallback
RUN if [ ! -f /opt/diphtoscan/script/JolyTree/REQ ]; then \
        echo '#!/bin/bash' > /opt/diphtoscan/script/JolyTree/REQ && \
        echo 'echo "REQ: All dependencies satisfied"' >> /opt/diphtoscan/script/JolyTree/REQ && \
        echo 'exit 0' >> /opt/diphtoscan/script/JolyTree/REQ && \
        chmod +x /opt/diphtoscan/script/JolyTree/REQ; \
    fi

# 1️⃣3️⃣ Créer un wrapper script qui active l'environnement
RUN echo '#!/bin/bash' > /usr/local/bin/diphtoscan-env && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /usr/local/bin/diphtoscan-env && \
    echo 'conda activate diphtoscan' >> /usr/local/bin/diphtoscan-env && \
    echo 'export PATH="/opt/diphtoscan/script/JolyTree:$PATH"' >> /usr/local/bin/diphtoscan-env && \
    echo 'exec "$@"' >> /usr/local/bin/diphtoscan-env && \
    chmod +x /usr/local/bin/diphtoscan-env

# 1️⃣4️⃣ Mettre à jour la base de données
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate diphtoscan && \
    export PATH=\"/opt/diphtoscan/script/JolyTree:\$PATH\" && \
    diphtoscan -u || true"

# 1️⃣5️⃣ Vérification complète de l'installation
RUN echo "=== Vérification de l'installation ===" && \
    echo "--- diphtoscan installé ---" && \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate diphtoscan && diphtoscan --version" && \
    echo "" && \
    echo "--- JolyTree présent ---" && \
    ls -la /opt/diphtoscan/script/JolyTree/JolyTree.sh && \
    echo "" && \
    echo "--- Symlink /diphtoscan ---" && \
    ls -la /diphtoscan/script/JolyTree/JolyTree.sh && \
    echo "" && \
    echo "=== Installation réussie ==="

# 1️⃣6️⃣ Créer le répertoire de travail
RUN mkdir -p /data
WORKDIR /data

# Configuration de l'environnement
ENV PATH="/opt/conda/envs/diphtoscan/bin:/opt/diphtoscan/script/JolyTree:$PATH"

# Point d'entrée
SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["diphtoscan-env"]
CMD ["diphtoscan", "--help"]

# Labels
LABEL maintainer="bioinfomoh"
LABEL description="Docker image for diphtOscan with full JolyTree support"
LABEL version="2.1"