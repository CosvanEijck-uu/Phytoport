#!/bin/bash
set -euo pipefail

source "$(dirname "$0")/input_utils.sh"

echo "===================================================="
echo " _____  _           _        _____           _  "
echo "|  __ \| |         | |      |  __ \         | |  "
echo "| |__) | |__  _   _| |_ ___ | |__) |__  _ __| |_ "
echo "|  ___/| '_ \| | | | __/ _ \|  ___/ _ \| '__| __|"
echo "| |    | | | | |_| | || (_) | |  | (_) | |  | |_ "
echo "|_|    |_| |_|\__, |\__\___/|_|   \___/|_|   \__|"
echo "                __/ |                            "
echo "              |___/                              "
echo "===================================================="
echo
echo "Welcome to PhytoPort — your gateway to projecting Arabidopsis knowledge to other plants."
echo

# -----------------------------
# Step 1: Second organism UPID
# -----------------------------
step1_get_second_organism_upid() {
    while true; do
        prompt_input "Please enter the UniProt proteome UPID of the comparison organism (e.g., UP000005640): " organism2_upid false
        if [[ "$organism2_upid" =~ ^UP[0-9]{6,9}$ ]]; then
            break
        fi
        echo "Invalid UPID format. It should look like: UP followed by 6–9 digits, e.g., UP000005640"
    done
    echo "Using second organism UPID: $organism2_upid"
    echo
}

# -----------------------------
# Step 2: GEO identifier (skippable)
# -----------------------------
step2_get_geo_id() {
    echo
    echo "============================================================"
    echo "📘 GEO Dataset Prerequisites"
    echo "============================================================"
    echo "Before continuing, please ensure you have a valid GEO dataset ID."
    echo
    echo "🔹 Accepted GEO Identifier Formats:"
    echo "   - GSE##### → Series (e.g., an experiment or project)"
    echo "   - GSM##### → Sample (e.g., an individual dataset within a series)"
    echo
    echo "🔹 Requirements for the dataset:"
    echo "   1. The dataset must be publicly available on NCBI GEO."
    echo "   2. It should correspond to gene expression or single-cell data."
    echo "   3. It should be in the 10X Genomics format:"
    echo "      <10x_directory>/ "
    echo "      ├── matrix.mtx[.gz] "
    echo "      ├── features.tsv[.gz]  or  genes.tsv[.gz] "
    echo "      └── barcodes.tsv[.gz] "
    echo "   4. IDs must map to EnsemblPlants GeneIDs via UniProt."
    echo "============================================================"
    echo
    echo "💡 Tip: Enter 'skip' to skip downloading or using a GEO dataset."

    while true; do
        prompt_input "Please enter a GEO identifier (e.g., GSE12345) or 'skip': " geo_id false
        if [[ "$geo_id" == "skip" ]]; then
            echo "⚠️  Skipping GEO dataset step."
            geo_id=""
            break
        elif [[ "$geo_id" =~ ^GS[EM][0-9]+$ ]]; then
            echo "✅ Using GEO identifier: $geo_id"
            break
        else
            echo "❌ Invalid GEO identifier format."
            echo "   It should look like: GSE12345 (series) or GSM67890 (sample), or 'skip'"
        fi
    done
    echo
}
# -----------------------------
# Step 3: Target proteins
# -----------------------------
step3_get_target_proteins() {
    while true; do
        read -p "Enter target proteins separated by spaces (e.g., HY5 SPA* COP1): " -a target_proteins
        if [ ${#target_proteins[@]} -gt 0 ]; then
            break
        fi
        echo "You must enter at least one target protein."
    done

    for i in "${!target_proteins[@]}"; do
        target_proteins[i]=$(echo "${target_proteins[$i]}" | tr '[:lower:]' '[:upper:]')
    done
}

# -----------------------------
# Step 4: Write config.yaml
# -----------------------------
step4_write_config() {
    config_file="config.yaml"

    to_yaml_list() {
        local arr=("$@")
        local result="["
        for el in "${arr[@]}"; do
            el=$(echo "$el" | sed 's/"/\\"/g')
            result+="\"$el\", "
        done
        result="${result%, }]"
        echo "$result"
    }

    inline_targets=$(to_yaml_list "${target_proteins[@]}")

    cat > "$config_file" << EOF
organism2_upid: "$organism2_upid"
geo_identifier: "$geo_id"
target_proteins: $inline_targets
structure_proteins: $inline_targets
EOF

    echo
    echo "Config written to $config_file:"
    cat "$config_file"
    echo
}

# -----------------------------
# Step 5: Run Snakemake
# -----------------------------
step5_run_snakemake() {
    echo "Running Snakemake workflow..."
    exec snakemake --snakefile /workspace/Snakefile --cores all --printshellcmds --latency-wait 20 --rerun-incomplete
}

# -----------------------------
# Main
# -----------------------------
start_pipeline() {
    step1_get_second_organism_upid
    step2_get_geo_id
    step3_get_target_proteins
    step4_write_config
    step5_run_snakemake
}

start_pipeline
