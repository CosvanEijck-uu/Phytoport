#!/usr/bin/env bash
set -euo pipefail

### ---------------------------------------------------------
###  1. Resolve absolute host paths (critical for Apptainer)
### ---------------------------------------------------------
# Results directory (same as Docker bind)
HOST_RESULTS="$(readlink -f ../Results)"
mkdir -p "$HOST_RESULTS"

# InterProScan extracted folder on the host
HOST_CACHE_DIR="$(pwd)/.interproscan_cache"
IPR_ACTUAL_HOST_DIR_NAME="interproscan-5.75-106.0"
HOST_IPR_DATA="${HOST_CACHE_DIR}/${IPR_ACTUAL_HOST_DIR_NAME}"

if [ ! -d "$HOST_IPR_DATA" ]; then
    echo "❌ ERROR: InterProScan data not found at $HOST_IPR_DATA"
    echo "Please extract it first."
    exit 1
fi

# Container paths (fixed, matches your Snakemake scripts)
CONTAINER_IPR_PATH="/opt/interproscan/interproscan-5.75-106.0-64-bit"
CONTAINER_OUTPUT_PATH="/workspace/output"

### ---------------------------------------------------------
###  2. Run Apptainer with Docker-equivalent behavior
### ---------------------------------------------------------
apptainer exec \
    --bind "${HOST_RESULTS}:${CONTAINER_OUTPUT_PATH}" \
    --bind "${HOST_IPR_DATA}:${CONTAINER_IPR_PATH}" \
    --pwd /workspace \
    --writable-tmpfs \
    phytoport.sif \
    bash run_boltz2.sh
