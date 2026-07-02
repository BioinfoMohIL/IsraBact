# R image for the --lang r path: ggtree + ggplot2 + ape.
# Build from the repo root:
#   docker build -f Dockerfile.R -t bioinfomoh/amr-heatmap-r:latest .
# (Bioconductor base already ships BiocManager; ggtree is a Bioconductor pkg.)
FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN R -e 'install.packages(c("ggplot2","ape","aplot"), repos="https://cloud.r-project.org")' \
 && R -e 'BiocManager::install("ggtree", ask=FALSE, update=FALSE)'

COPY scripts/ /opt/amr-heatmap/scripts/
RUN chmod +x /opt/amr-heatmap/scripts/*.R
ENV PATH=/opt/amr-heatmap/scripts:$PATH

RUN R -e 'library(ggtree); library(ggplot2); library(ape); library(aplot); cat("R deps ok\n")'

WORKDIR /data