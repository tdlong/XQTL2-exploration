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

# 1) Adaptive list-format estimation (fresh run)
Rscript scripts/production/run_haplotype_estimation_list_format.R \
  "${CHR}" adaptive 4 "${OUTPUT_DIR}" "${PARAM_FILE}"

# 2) Smooth h4 from that adaptive list-format (writes both adapt/smooth reshaped files)
Rscript scripts/production/create_smooth_haplotype_estimator_list_format.R \
  "${CHR}" "${PARAM_FILE}" "${OUTPUT_DIR}"

echo "[INFO] $(date) done ${CHR}"
