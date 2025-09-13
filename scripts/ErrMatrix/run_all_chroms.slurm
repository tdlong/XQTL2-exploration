#!/usr/bin/env bash

#SBATCH --job-name=base_var_wide_all
#SBATCH --output=logs/base_var_wide_all_%A_%a.out
#SBATCH --error=logs/base_var_wide_all_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-5

set -euo pipefail

# Chromosomes to process
chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
CHR=${chrs[$((SLURM_ARRAY_TASK_ID-1))]}

# Use JUICE parameters (not ZINC2) for ErrMatrix testing
PARAM_FILE=helpfiles/JUICE_haplotype_parameters.R
OUTPUT_DIR=process/JUICE

mkdir -p logs

echo "[INFO] $(date) starting ${CHR} with BASE_VAR_WIDE.R"
module load R/4.4.2

# Run BASE_VAR_WIDE.R in production mode (non-debug, non-verbose)
# This will run the complete workflow: adaptive estimation + smoothing
Rscript scripts/ErrMatrix/BASE_VAR_WIDE.R \
  "${CHR}" adaptive 4 "${OUTPUT_DIR}" "${PARAM_FILE}" --nonverbose

echo "[INFO] $(date) done ${CHR}"
