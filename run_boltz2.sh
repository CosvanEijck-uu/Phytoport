#!/bin/bash
set -euo pipefail

echo
echo "===================================================="
echo "     ⚛️  Boltz-2 Structure Prediction Stage"
echo "===================================================="
echo
echo "This step performs deep-learning-based protein structure predictions."
echo "⚠️  It is computationally intensive and may take a long time, even on GPU-enabled systems."
echo

config_file="config.yaml"

if [ ! -f "$config_file" ]; then
    echo "❌ Error: $config_file not found. Please run the main PhytoPort pipeline first."
    exit 1
fi

# Load target proteins from config
target_proteins=$(yq '.target_proteins | join(" ")' "$config_file")

echo "Detected target proteins from main pipeline:"
echo "$target_proteins"
echo

# Let user choose structure proteins (default = same as targets)
read -p "Would you like to restrict structure prediction to a subset? [y/N]: " subset
if [[ "$subset" =~ ^[Yy]$ ]]; then
    read -p "Enter subset of proteins separated by spaces: " -a subset_proteins
    echo "Using subset: ${subset_proteins[*]}"

    # Convert bash array to YAML list string
    inline_subset="["
    for el in "${subset_proteins[@]}"; do
        inline_subset+="\"$el\", "
    done
    inline_subset="${inline_subset%, }]"  # Remove trailing comma

    # Update YAML in-place with yq v4+
    yq -i -y ".structure_proteins = $inline_subset" "$config_file"
else
    # Copy target_proteins to structure_proteins
    yq -i -y ".structure_proteins = .target_proteins" "$config_file"
fi

echo
echo "✅ Updated structure_proteins in config.yaml:"
yq -y '.structure_proteins' "$config_file"
echo

echo "Starting Boltz-2 Snakemake workflow..."
snakemake --snakefile ./Snakefile --cores all --printshellcmds --latency-wait 20 boltz_predict
