#!/bin/bash

# --- Configuration ---
IMAGE_NAME="cvaneijck/phytoport:latest"
SIF_FILE="phytoport.sif"

# SIF Build Temporary Directory (for apptainer build process)
TMP_DIR="${PWD}/.tmp_apptainer_build_$$"

# InterProScan Data Cache Directory (for the extracted 39GB data)
HOST_CACHE_DIR="${PWD}/.interproscan_cache"

# Large temporary working directory for extraction (using /tmp)
HOST_WORK_DIR="/tmp/.apptainer_work_dir_$$"

# Internal paths defined in the container's entrypoint
IPR_ARCHIVE_PATH="/opt/interproscan/interproscan-5.75-106.0-64-bit.tar.gz"
IPR_TARGET_DIR="/opt/interproscan/interproscan-5.75-106.0-64-bit"

# ACTUAL extracted folder name from tar output (observed as 'interproscan-5.75-106.0')
IPR_ACTUAL_HOST_DIR_NAME="interproscan-5.75-106.0"
# --- End Configuration ---

echo "Starting Apptainer build process for image: ${IMAGE_NAME}"
echo "Output SIF file: ${SIF_FILE}"
echo "Extracted data will be stored in: ${HOST_CACHE_DIR}"

# --- Part 1: SIF Build (Pulling the Docker Image) ---

# Check if Apptainer is installed
if ! command -v apptainer &> /dev/null; then
    echo "🚨 ERROR: Apptainer command not found."
    echo "Please install Apptainer to run this script."
    exit 1
fi

# 1. Create the temporary directory for the SIF build (using absolute path for robustness)
echo "Creating temporary build directory: ${TMP_DIR}"
mkdir -p "${TMP_DIR}"

# 2. Set TMPDIR environment variable for SIF building
export TMPDIR="${TMP_DIR}"
echo "Setting TMPDIR to: ${TMPDIR}"

# 3. Pull the Docker image and create the SIF file
echo "Pulling Docker image and building SIF file. This may take a while..."
apptainer pull "${SIF_FILE}" "docker://${IMAGE_NAME}"

# Check the exit status of the pull command
if [ $? -ne 0 ]; then
    echo "❌ ERROR: Apptainer pull failed. Check the error messages above."
    unset TMPDIR
    rm -rf "${TMP_DIR}"
    exit 1
fi

# 4. Cleanup SIF Build Temp Dir
echo "Unsetting TMPDIR and cleaning up SIF build temp directory: ${TMP_DIR}"
unset TMPDIR
rm -rf "${TMP_DIR}"

# --- Part 2: Data Extraction (One-Time Setup) ---

echo "--- Starting One-Time Data Extraction ---"

# Define the full path to the expected extracted folder
IPR_HOST_DIR="${HOST_CACHE_DIR}/${IPR_ACTUAL_HOST_DIR_NAME}"

# 5. Check if the extracted directory already exists
if [ -d "${IPR_HOST_DIR}" ]; then
    echo "✅ InterProScan data already found at ${IPR_HOST_DIR}. Skipping extraction."
else
    echo "⚠️ InterProScan data not found. Proceeding with extraction (approx 39GB)."

    # Create host cache directory
    echo "Creating host data cache directory: ${HOST_CACHE_DIR}"
    mkdir -p "${HOST_CACHE_DIR}"

    # Create host work directory in /tmp
    echo "Creating host work directory for extraction buffers: ${HOST_WORK_DIR}"
    mkdir -p "${HOST_WORK_DIR}"

    # 6. Extract the archive to the host directory using --workdir in /tmp
    apptainer exec \
        --bind "${HOST_CACHE_DIR}":/tmp/writable_extract \
        --workdir "${HOST_WORK_DIR}" \
        "${SIF_FILE}" \
        /bin/bash -c "tar -pxvzf ${IPR_ARCHIVE_PATH} -C /tmp/writable_extract"

    # Check the exit status of the extraction command
    if [ $? -ne 0 ]; then
        echo "❌ ERROR: Data extraction failed. Check disk space on /tmp or your main partition."
        rm -rf "${HOST_WORK_DIR}"
        exit 1
    fi

    # Clean up the temporary work directory created in /tmp
    echo "Cleaning up temporary work directory: ${HOST_WORK_DIR}"
    rm -rf "${HOST_WORK_DIR}"
fi

# 7. Final Cleanup and Summary

# Ensure the SIF file is owned by the current user
echo "Setting ownership of ${SIF_FILE} to current user."
chown $(id -u):$(id -g) "${SIF_FILE}"

# Define the final usage paths for the user
CONTAINER_IPR_PATH="${IPR_TARGET_DIR}"
# Define the final usage paths for the user
# ... (IPR_HOST_DIR and CONTAINER_IPR_PATH remain the same) ...
CONTAINER_OUTPUT_PATH="/workspace/output"

echo "✅ Installation and data setup complete! The pipeline is ready to run."
echo "--------------------------------------------------------"
echo "To run your pipeline, use the following command:"
echo "       bash run_phytoport.sh"
echo "--------------------------------------------------------"
