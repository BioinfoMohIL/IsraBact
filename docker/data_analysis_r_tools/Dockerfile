# Base image with R and tidyverse pre-installed
FROM rocker/tidyverse:latest

# Install additional system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('pheatmap', 'patchwork', 'circlize', 'janitor', 'RColorBrewer', 'argparse'))"
