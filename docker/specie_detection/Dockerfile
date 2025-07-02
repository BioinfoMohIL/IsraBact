FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required dependencies
RUN apt-get update && apt-get install -y \
    wget \
    tar \
    g++ \
    make \
    rsync \
    perl \
    build-essential \
    libz-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Download and extract Kraken2
RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.14.tar.gz -O kraken2.tar.gz && \
    tar -xzf kraken2.tar.gz && \
    mv kraken2-2.14 kraken2 && \
    rm kraken2.tar.gz

# Install Kraken2
WORKDIR /opt/kraken2
RUN chmod +x install_kraken2.sh && \
    ./install_kraken2.sh /usr/local/bin

    
# Add Kraken2 to PATH
ENV PATH="/usr/local/bin:${PATH}"
    
    
WORKDIR /app  
COPY . /app/

RUN chmod +x /app/specie_detection.sh

CMD ["bash"]

