# Use an official lightweight Python base image
FROM python:3.10-slim

# Metadata
LABEL maintainer="bioinfomoh"
LABEL version="1.0"
LABEL description="Salmonella/EColi Virulence Finder Tool"

# Install required system packages (BLAST+, gzip, etc.)
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set work directory
WORKDIR /app

# Copy script and data folder into container
COPY virulence_finder.py /app/virulence_finder.py
COPY data/ /app/data/

# Install Python dependencies
RUN pip install --no-cache-dir biopython pandas jinja2

# Add a simple CLI entry point
RUN echo '#!/bin/bash\npython /app/virulence_finder.py "$@"' > /bin/virulence_finder \
    && chmod +x /bin/virulence_finder
