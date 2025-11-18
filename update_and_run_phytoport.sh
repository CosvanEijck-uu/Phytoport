#!/bin/bash

set -euo pipefail

# -----------------------------
# Check if Docker is running
# -----------------------------
if ! docker info >/dev/null 2>&1; then
  echo "Docker does not seem to be running. Please start Docker and try again."
  exit 1
fi

# -----------------------------
# Build main pipeline Docker image
# -----------------------------
echo "Building Docker image 'phytoport' from Dockerfile in $(pwd)..."
docker build -t phytoport:latest .

# -----------------------------
# Run the main pipeline container
# -----------------------------
# echo "Running 'phytoport' container..."
# docker run -it --rm --shm-size=4g -v "$(pwd)/../Results/:/workspace/output" phytoport_snakemake:latest

