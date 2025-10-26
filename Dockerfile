# Use a lightweight base image with Ubuntu
FROM ubuntu:22.04

# Metadata
LABEL maintainer="Your Name <you@example.com>" \
      description="Docker image for downloading and compressing RNA-seq FASTQ files using SRA Toolkit (fasterq-dump)."

# Install dependencies
RUN apt-get update && apt-get install -y \
    sra-toolkit \
    pigz \
    wget \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /data

# Environment variables
ENV SRA_TOOLKIT_HOME=/usr/bin
ENV PATH=$PATH:$SRA_TOOLKIT_HOME

# Default command — can be overridden
ENTRYPOINT ["/bin/bash", "-c"]
