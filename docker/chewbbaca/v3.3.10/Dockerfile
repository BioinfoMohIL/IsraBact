FROM ummidock/chewbbaca:v3.3.10

RUN apt-get update && \
    apt-get install -y \
        unzip \
        zip \
        tree \
        jq \
        rsync \
        curl \
        wget \
        apt-transport-https \
        ca-certificates \
        gnupg && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" > \
        /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && \
    apt-get install -y google-cloud-cli && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*