FROM mcr.microsoft.com/devcontainers/python:3.12

# Runtime prerequisites for Nextflow and this pipeline
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow

# Python dependencies used by scripts in this repository
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
      torch \
      torchvision \
      celltypist \
      leidenalg \
      scvi-tools \
      scanpy \
      scikit-learn

# Optional dependency: keep build resilient if RAPIDS wheels are unavailable
RUN pip install --no-cache-dir rapids-singlecell || true
