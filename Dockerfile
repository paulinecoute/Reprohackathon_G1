# Use a lightweight base image with Ubuntu
FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    sra-toolkit \
    pigz \
    wget \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Environment variables
ENV SRA_TOOLKIT_HOME=/usr/bin
ENV PATH=$PATH:$SRA_TOOLKIT_HOME

# Default command — can be overridden
ENTRYPOINT ["/bin/bash", "-c"]
