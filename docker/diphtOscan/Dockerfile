FROM ubuntu:22.04 
ENV DEBIAN_FRONTEND=noninteractive 

# 1️⃣ Dépendances système 
RUN apt-get update && apt-get install -y \ 
    wget curl git bzip2 build-essential unzip \ 
    hmmer ncbi-blast+ mash \ 
    && apt-get clean && rm -rf /var/lib/apt/lists/* 

# 2️⃣ Installer Miniconda 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \ 
    && bash /tmp/miniconda.sh -b -p /opt/conda \ && rm /tmp/miniconda.sh ENV PATH="/opt/conda/bin:$PATH" 

# 3️⃣ Cloner diphtoscan 
RUN git clone https://gitlab.pasteur.fr/BEBP/diphtoscan.git /opt/diphtoscan WORKDIR /opt/diphtoscan 

#4️⃣ Accepter les Terms of Service 
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \ 
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 

# 5️⃣ Créer l'environnement conda 
RUN conda env create -f environment.yml 

# 6️⃣ Installer les dépendances manquantes 
RUN /bin/bash -c "source activate diphtoscan && \ 
    conda install -y -c bioconda ncbi-amrfinderplus infernal prodigal" 
    
# 7️⃣ Installer diphtoscan dans l'environnement 
RUN /bin/bash -c "source activate diphtoscan && python -m pip install . --no-deps" 

# 8️⃣ Ajouter les binaires au PATH global 
RUN ln -s /opt/conda/envs/diphtoscan/bin/diphtoscan /usr/local/bin/diphtoscan \ 
    && ln -s /opt/conda/envs/diphtoscan/bin/amrfinder /usr/local/bin/amrfinder \ 
    && ln -s /opt/conda/envs/diphtoscan/bin/cmsearch /usr/local/bin/cmsearch \ 
    && ln -s /opt/conda/envs/diphtoscan/bin/prodigal /usr/local/bin/prodigal 

# 9️⃣ Mettre à jour la base 
RUN /bin/bash -c "source activate diphtoscan && diphtoscan -u"