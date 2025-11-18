# -----------------------------
# Base image
# -----------------------------
FROM ubuntu:22.04

# -----------------------------
# Environment variables
# -----------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH
ENV DISPLAY=:99

# -----------------------------
# Install system dependencies
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git bzip2 unzip tar gzip \
    build-essential gfortran cmake make \
    python3 python3-pip python3-venv python3-distutils \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libbz2-dev liblzma-dev zlib1g-dev \
    libx11-dev libxt-dev libpng-dev libjpeg-dev \
    libreadline-dev libncurses5-dev libncursesw5-dev \
    tzdata libpcre2-dev \
    perl openjdk-11-jdk gnupg ca-certificates \
    sudo vim less tree \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Install Miniconda
# -----------------------------
RUN wget --progress=dot:giga https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean --all -y

# -----------------------------
# Accept Conda Terms
# -----------------------------
RUN /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# -----------------------------
# Create bioinformatics environment (Python 3.11 + R 4.4.3)
# -----------------------------
RUN /opt/conda/bin/conda create -n bioenv python=3.11 r-base=4.4.3 r-essentials -c conda-forge -y

SHELL ["/opt/conda/bin/conda", "run", "-n", "bioenv", "/bin/bash", "-c"]

# -----------------------------
# Install Python packages
# -----------------------------
COPY requirements.txt /workspace/requirements.txt
WORKDIR /workspace
RUN python -m pip install --upgrade pip && \
    python -m pip install -r requirements.txt

# Install Snakemake inside bioenv
RUN conda install -c bioconda -c conda-forge snakemake -y

SHELL ["/bin/bash", "-c"]

# -----------------------------
# ETE3 environment (Python 3.6 + PyQt4)
# -----------------------------
RUN /opt/conda/bin/conda create -n ete3_env python=3.6 -y
SHELL ["/opt/conda/bin/conda", "run", "-n", "ete3_env", "/bin/bash", "-c"]
RUN conda install -c bioconda ete3 -y
SHELL ["/bin/bash", "-c"]

# Add ete3_env to PATH so Snakemake can access it
ENV PATH=/opt/conda/envs/ete3_env/bin:$PATH

# -----------------------------
# InterProScan
# -----------------------------
WORKDIR /opt/interproscan
RUN wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz && \
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz.md5 && \
    md5sum -c interproscan-5.75-106.0-64-bit.tar.gz.md5

ENV INTERPROSCAN_PATH=/opt/interproscan/interproscan-5.75-106.0
RUN /opt/conda/bin/conda run -n bioenv bash -c "conda env config vars set JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64"

# -----------------------------
# OrthoFinder
# -----------------------------
RUN /opt/conda/bin/conda create -n of3_env python=3.11 -y && \
    /opt/conda/bin/conda run -n of3_env pip install git+https://github.com/OrthoFinder/OrthoFinder.git

# -----------------------------
# Add conda & OrthoFinder envs to PATH
# -----------------------------
ENV PATH=/opt/conda/envs/bioenv/bin:/opt/of3_env/bin:$PATH

# -----------------------------
# MusiteDeep (Python 2.7)
# -----------------------------
RUN /opt/conda/bin/conda create -n musitedeep_env python=2.7 -y && \
    /opt/conda/bin/conda run -n musitedeep_env pip install --upgrade pip && \
    /opt/conda/bin/conda run -n musitedeep_env pip install numpy==1.16.6 scipy==1.2.3 pandas==0.24.2 h5py==2.10.0 keras==2.1.2 tensorflow==1.3.0

RUN git clone https://github.com/duolinwang/MusiteDeep.git /MusiteDeep_Keras2.0
WORKDIR /MusiteDeep_Keras2.0/MusiteDeep_Keras2.0/MusiteDeep

# -----------------------------
# Boltz-2 (Python 3.9, CUDA support)
# -----------------------------
RUN /opt/conda/bin/conda create -n boltz2_env python=3.9 -y && \
    /opt/conda/bin/conda run -n boltz2_env pip install --upgrade pip && \
    /opt/conda/bin/conda run -n boltz2_env pip install "boltz[cuda]" -U

# Optionally verify installation
RUN /opt/conda/bin/conda run -n boltz2_env python -m boltz --help || true

# -----------------------------
# Install R packages
# -----------------------------
COPY installed_packages.csv /tmp/installed_packages.csv
RUN conda install -n bioenv -c conda-forge libpng zlib -y
RUN /opt/conda/bin/conda run -n bioenv Rscript -e '\
  pkgs <- read.csv("/tmp/installed_packages.csv", stringsAsFactors = FALSE)$Package; \
  repos <- c(CRAN="https://cloud.r-project.org"); \
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos=repos); \
  BiocManager::install(pkgs, ask=FALSE, update=TRUE, Ncpus=parallel::detectCores()); \
  message("✅ Installed R packages: ", length(pkgs)) \
  '

# -----------------------------
# Copy pipeline scripts
# -----------------------------
COPY . /workspace
WORKDIR /workspace

# -----------------------------
# Entrypoint for Snakemake
# -----------------------------
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
