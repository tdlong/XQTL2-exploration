#!/usr/bin/env bash

#SBATCH --job-name=haps_list_all
#SBATCH --output=logs/haps_list_all_%A_%a.out
#SBATCH --error=logs/haps_list_all_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-5

set -euo pipefail

# Chromosomes to process
chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
CHR=${chrs[$((SLURM_ARRAY_TASK_ID-1))]}

PARAM_FILE=helpfiles/ZINC2_haplotype_parameters.R
OUTPUT_DIR=process/ZINC2

mkdir -p logs

echo "[INFO] $(date) starting ${CHR}"
module load R/4.4.2

# Run complete workflow (adaptive estimation + smoothing) in one call
# Using non-verbose mode for production
Rscript scripts/production/complete_haplotype_workflow.R \
  "${CHR}" adaptive 4 "${OUTPUT_DIR}" "${PARAM_FILE}" --nonverbose

echo "[INFO] $(date) done ${CHR}"
