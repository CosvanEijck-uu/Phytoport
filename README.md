# PhytoPort Pipeline Documentation 🧬

PhytoPort is a bioinformatics pipeline for **cross-species protein analysis**. It facilitates the extraction of orthologs, comparison of protein domains and phosphorylation sites, gene expression analysis, Protein-Protein Interaction (PPI) studies, and management of reference proteomes.

-----

## 🛠️ Pipeline Execution & Environment

The pipeline is primarily managed using **Snakemake** and distributed via **Apptainer** or **local installation**.

  * `README.md`: General overview and main entry point documentation.
  * `phytoport.sif`: The Singularity Image Format file for container execution.

### Shell Scripts

| Script | Description |
| :--- | :--- |
| `initialyse.sh` | **Main wrapper script** that provides a CLI. It prompts the user for the **comparison organism** (by name/UPID), **target proteins**, and **structure prediction targets**. It generates `config.yaml` and executes the pipeline via the `Snakefile`. |
| `entrypoint.sh` | Configures the **Docker entry point**. Handles container startup, volume mounting, and initiates the pipeline execution. |
| `build_phytoport.sh` | This script builds an Apptainer SIF image and extracts required InterProScan data to the host system. |
| `run_boltz2.sh` | Sets up and runs **Boltz-2** deep-learning protein structure predictions. Allows users to optionally restrict targets and updates `config.yaml` before launching the Snakemake workflow. |
| `input_utils.sh` | Utility functions for parsing, validating, and preparing input for the pipeline. |

-----

## 💻 Python Modules

The pipeline's core functionality is organized into modular Python packages: `fetch`, `extract`, `compare`, `summarize`, and `visualize`.

### 📂 `fetch` (Data Retrieval)

| Script | Functionality |
| :--- | :--- |
| `proteomes.py` | Downloads **reference proteomes** from UniProt using organism name or UPID, supporting multiple organisms. |
| `proteins.py` | Fetches **FASTA files** for specified proteins from UniProt, preparing them for PPI analysis (e.g., creating pairwise FASTA files). |
| `geo_dataset.py` | Downloads the first **10X-style GEO files** for a given GSE or GSM ID, using FTP/HTTPS. |

### 🔍 `extract` (Ortholog and Tree Processing)

| Script | Functionality |
| :--- | :--- |
| `target_orthologs.py` | Extracts target orthologs from an **OrthoFinder TSV**. Matches genes by target symbols (supports wildcards like `SPA*`), organizes results, and outputs filtered TSV files. |
| `orthogroup_trees.py` | Extracts specific **orthogroup trees** from a Newick file based on a list of Orthogroup IDs. |

### ⚖️ `compare` (Feature Comparison)

| Script | Functionality |
| :--- | :--- |
| `domains.py` | Compares predicted **domains/motifs** between two protein sequences from a TSV, identifying common and unique features. |
| `MSA_ps_sites.py` | Compares predicted **phosphorylation (PS) sites** between two protein sequences, identifying common/unique sites by MSA column, with optional score filtering. |
| `structures.py` | This script submits two CIF/mmCIF files to the RCSB FATCAT alignment API and prints a summary of the structural alignment. |

### 📊 `summarize` (Prediction Processing)

| Script | Functionality |
| :--- | :--- |
| `interproscan.py` | Summarizes raw **InterProScan** output for use in downstream analysis (e.g., domain comparison). |
| `musitedeep.py` | Summarizes **MusiteDeep** phosphorylation site predictions, including filtering by score cutoffs. |

### 🖼️ `visualize` (Results Presentation)

| Script | Functionality |
| :--- | :--- |
| `ps_sites.py` | **Visualizes PS site scores** on protein sequences from a FASTA file using color and font size to represent per-residue score values. |
| `tree.py` | Renders a **Newick phylogenetic tree** with optional sequence and domain visualization. Automatically caps branch lengths to mitigate outlier skew. |

-----

## 🧪 Testing

  * The pipeline has been validated through a **full case study execution**.

-----

## 📚 External Dependencies & Resources

  * **OrthoFinder:** [Installation Guide](https://github.com/OrthoFinder/OrthoFinder?tab=readme-ov-file#installation)
    - Reference: Emms, D. M., Liu, Y., Belcher, L., Holmes, J., & Kelly, S. (2025). OrthoFinder: Scalable phylogenetic orthology inference for comparative genomics (p. 2025.07.15.664860). bioRxiv. https://doi.org/10.1101/2025.07.15.664860
  * **InterProScan:** [Documentation](https://interpro-documentation.readthedocs.io/en/latest/interproscan.html)
    - Reference: Jones, P., Binns, D., Chang, H.-Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G., Pesseat, S., Quinn, A. F., Sangrador-Vegas, A., Scheremetjew, M., Yong, S.-Y., Lopez, R., & Hunter, S. (2014). InterProScan 5: Genome-scale protein function classification. Bioinformatics, 30(9), 1236–1240. https://doi.org/10.1093/bioinformatics/btu031
  * **MusiteDeep:** [GitHub Repository](https://github.com/duolinwang/MusiteDeep)
    - Reference: Wang, D., Zeng, S., Xu, C., Qiu, W., Liang, Y., Joshi, T., & Xu, D. (2017). MusiteDeep: A deep-learning framework for general and kinase-specific phosphorylation site prediction. Bioinformatics (Oxford, England), 33(24), 3909–3916. https://doi.org/10.1093/bioinformatics/btx496
  * **Boltz-2:** [GitHub Repository](https://github.com/jwohlwend/boltz/tree/main)
    - Reference: Passaro, S., Corso, G., Wohlwend, J., Reveiz, M., Thaler, S., Somnath, V. R., Getz, N., Portnoi, T., Roy, J., Stark, H., Kwabi-Addo, D., Beaini, D., Jaakkola, T., & Barzilay, R. (2025). Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction (p. 2025.06.14.659707). bioRxiv. https://doi.org/10.1101/2025.06.14.659707
  * **Seurat:** [Documentation](https://satijalab.org/seurat/)
    - Reference: Gribov, A., Sill, M., Lück, S., Rücker, F., Döhner, K., Bullinger, L., Benner, A., & Unwin, A. (2010). SEURAT: Visual analytics for the integrated analysis of microarray data. BMC Medical Genomics, 3(1), 21. https://doi.org/10.1186/1755-8794-3-21

-----

## 🚀 Getting Started (Manual Steps)

1.  Ensure you are on a **Linux** or **WSL** command line.

2.  Ensure **Apptainer** is installed and running.

3.  Run the installation script to build the image:

    ```bash
    bash build_phytoport.sh
    ```

4.  Run the execution script to run the image (This mounts the correct directories 🙂):

```bash
    bash run_phytoport.sh
```


