# 🔧 **Local Installation Guide**

This guide installs all components **directly on Ubuntu 22.04 or similar Linux**.

-----

## 1\. 📦 Install System Dependencies

```bash
sudo apt update
sudo apt install -y \
    wget curl git bzip2 unzip tar gzip \
    build-essential gfortran cmake make \
    python3 python3-pip python3-venv python3-distutils \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libbz2-dev liblzma-dev zlib1g-dev \
    libx11-dev libxt-dev libpng-dev libjpeg-dev \
    libreadline-dev libncurses5-dev libncursesw5-dev \
    tzdata libpcre2-dev \
    perl openjdk-11-jdk gnupg ca-certificates \
    sudo vim less tree xvfb
```

-----

## 2\. 🟦 Install **Miniconda**

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
rm miniconda.sh

echo 'export PATH=$HOME/miniconda/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Accept Conda terms (same as Dockerfile):

```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

-----

## 3\. 🧬 Create **bioenv** (Python 3.11 + R 4.4.3)

```bash
conda create -n bioenv python=3.11 r-base=4.4.3 r-essentials -c conda-forge -y
conda activate bioenv
```

### Install Python packages

Copy your `requirements.txt` locally first, then:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### Install Snakemake

```bash
conda install -c bioconda -c conda-forge snakemake -y
```

-----

## 4\. 🌳 Create **ete3\_env** (Python 3.6 + ETE3)

```bash
conda create -n ete3_env python=3.6 -y
conda activate ete3_env
conda install -c bioconda ete3 -y
conda deactivate
```

-----

## 5\. 🧬 Create **OrthoFinder 3 environment**

```bash
conda create -n of3_env python=3.11 -y
conda activate of3_env
pip install git+https://github.com/OrthoFinder/OrthoFinder.git
conda deactivate
```

-----

## 6\. 🔬 Create **musitedeep\_env** (Python 2.7)

```bash
conda create -n musitedeep_env python=2.7 -y
conda activate musitedeep_env

pip install --upgrade pip
pip install numpy==1.16.6 scipy==1.2.3 pandas==0.24.2 \
            h5py==2.10.0 keras==2.1.2 tensorflow==1.3.0
```

Clone MusiteDeep:

```bash
git clone https://github.com/duolinwang/MusiteDeep.git ~/MusiteDeep_Keras2.0
conda deactivate
```

-----

## 7\. ⚡ Create **boltz2\_env** (Python 3.9 + CUDA support)

```bash
conda create -n boltz2_env python=3.9 -y
conda activate boltz2_env

pip install --upgrade pip
pip install "boltz[cuda]" -U
conda deactivate
```

Optional test:

```bash
conda run -n boltz2_env python -m boltz --help
```

-----

## 8\. 🧬 Install **InterProScan 5.75-106.0**

```bash
sudo mkdir -p /opt/interproscan
cd /opt/interproscan

# Download
sudo wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz
sudo wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.75-106.0-64-bit.tar.gz.md5
```

### Extract

```bash
# Extract the tarball
sudo tar -xvzf interproscan-5.75-106.0-64-bit.tar.gz
```

### Configure

Set env variable (add to `~/.bashrc`):

```bash
echo 'export INTERPROSCAN_PATH=/opt/interproscan/interproscan-5.75-106.0' >> ~/.bashrc
# Set JAVA_HOME for InterProScan within bioenv
conda run -n bioenv conda env config vars set JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
```

-----

## 9\. 🧬 Install **R packages** from `installed_packages.csv`

```bash
# Ensure libpng and zlib are available in bioenv for R packages
conda install -n bioenv -c conda-forge libpng zlib -y

# Run the R install script inside bioenv
conda run -n bioenv Rscript -e '
pkgs <- read.csv("installed_packages.csv", stringsAsFactors = FALSE)$Package;
repos <- c(CRAN="https://cloud.r-project.org");
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos=repos);
BiocManager::install(pkgs, ask=FALSE, update=TRUE, Ncpus=parallel::detectCores());
message("✅ Installed R packages: ", length(pkgs))
'
```

-----

## 10\. 📁 Pipeline Scripts

Place your workflow files inside a working directory:

```
/workspace/
    Snakefile
    config/
    scripts/
    ...
```

-----

## 11\. ▶️ Running Snakemake (like entrypoint.sh)

Activate **bioenv**, then run. The paths for `ete3_env` and `of3_env` are needed if your Snakemake workflow invokes tools from these environments.

```bash
source ~/.bashrc
conda activate bioenv

# Add ete3 and OrthoFinder to PATH for Snakemake to find them
export PATH=$HOME/miniconda/envs/ete3_env/bin:$HOME/miniconda/envs/of3_env/bin:$PATH

snakemake --cores all
```

-----

## ✅ Installation Complete

| Environment      | Purpose                                                      |
| ---------------- | ------------------------------------------------------------ |
| `bioenv`         | Snakemake + Python 3.11 + R 4.4.3 pipeline environment       |
| `ete3_env`       | ETE3 tree visualization (PyQt4)                              |
| `of3_env`        | OrthoFinder 3                                                |
| `musitedeep_env` | MusiteDeep (TensorFlow 1.3 + Keras 2.1.2)                    |
| `boltz2_env`     | Boltz (CUDA-enabled)                                         |
| InterProScan     | Functional annotation tool                                   |
