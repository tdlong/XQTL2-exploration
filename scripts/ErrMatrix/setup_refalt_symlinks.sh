#!/usr/bin/env bash

# Setup symlinks for RefAlt files to enable haplotype estimation
# Usage: ./setup_refalt_symlinks.sh <dataset> <output_dir>
# Example: ./setup_refalt_symlinks.sh JUICE process/JUICE_adapt_h10

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <dataset> <output_dir>"
    echo "Example: $0 JUICE process/JUICE_adapt_h10"
    echo "Example: $0 ZINC2 process/ZINC2_adapt_h10"
    exit 1
fi

DATASET="$1"
OUTPUT_DIR="$2"

# Source directory where RefAlt files are located
SOURCE_DIR="process/${DATASET}"

echo "=== SETTING UP REFALT SYMLINKS ==="
echo "Dataset: ${DATASET}"
echo "Source directory: ${SOURCE_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Create symlinks for all RefAlt files
for chr in chrX chr2L chr2R chr3L chr3R; do
    source_file="${SOURCE_DIR}/RefAlt.${chr}.txt"
    symlink_file="${OUTPUT_DIR}/RefAlt.${chr}.txt"
    
    if [ -f "${source_file}" ]; then
        if [ -L "${symlink_file}" ] || [ -f "${symlink_file}" ]; then
            echo "Removing existing: ${symlink_file}"
            rm -f "${symlink_file}"
        fi
        
        echo "Creating symlink: ${symlink_file} -> ${source_file}"
        ln -s "../${DATASET}/RefAlt.${chr}.txt" "${symlink_file}"
    else
        echo "WARNING: Source file not found: ${source_file}"
    fi
done

echo ""
echo "=== SYMLINKS CREATED ==="
echo "You can now run the haplotype estimation:"
echo "sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/${DATASET}_haplotype_parameters.R ${OUTPUT_DIR} 10"

