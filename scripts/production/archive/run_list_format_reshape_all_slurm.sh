#!/usr/bin/env bash

#SBATCH --job-name=reshape_list_format
#SBATCH --output=slurm_logs/reshape_list_format_%A_%a.out
#SBATCH --error=slurm_logs/reshape_list_format_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-5

set -euo pipefail

# Chromosomes to process (index with SLURM_ARRAY_TASK_ID)
chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
CHR=${chrs[$((SLURM_ARRAY_TASK_ID-1))]}

# Paths (relative, cluster-safe)
PARAM_FILE=helpfiles/ZINC2_haplotype_parameters.R
OUTPUT_DIR=process/ZINC2
SCRIPT=scripts/production/create_smooth_haplotype_estimator_list_format.R

# Ensure log directory exists
mkdir -p slurm_logs

echo "[INFO] Starting list-format reshape for ${CHR}"
echo "[INFO] Using PARAM_FILE=${PARAM_FILE} OUTPUT_DIR=${OUTPUT_DIR}"

# Run once per chromosome. This script writes both:
# - adapt_h4/R.haps.<chr>.out.rds
# - smooth_h4/R.haps.<chr>.out.rds
Rscript "${SCRIPT}" "${CHR}" "${PARAM_FILE}" "${OUTPUT_DIR}"

echo "[INFO] Completed ${CHR}"

