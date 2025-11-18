#!/bin/bash
set -euo pipefail

# -----------------------------
# Environment setup
# -----------------------------
# NOTE: The extraction logic is REMOVED because the data is now pre-extracted
# and mounted from the host machine during Apptainer execution.
# We rely on the Apptainer --bind flag to place the extracted data at:
# /opt/interproscan/interproscan-5.75-106.0-64-bit

# Add InterProScan and OrthoFinder to PATH
# CRITICAL: Ensure INTERPROSCAN_PATH points to the correct directory structure.
# If the path is relative to /opt/interproscan, you might need to set it here.
# Assuming the image sets INTERPROSCAN_PATH correctly, or that the path below is enough:
export PATH="/opt/interproscan/interproscan-5.75-106.0-64-bit:/opt/of3_env/bin:$PATH"

# Activate bioenv
source /opt/conda/bin/activate bioenv

# -----------------------------
# Default behavior
# -----------------------------
if [ $# -eq 0 ]; then
    echo "No arguments provided. Running PhytoPort wrapper script..."
    # Call the wrapper script
    exec ./initialyse.sh
else
    echo "Running command: $*"
    exec "$@"
fi
