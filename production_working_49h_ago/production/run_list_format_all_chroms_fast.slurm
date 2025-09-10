#!/usr/bin/env bash

#SBATCH --job-name=haps_list_fast
#SBATCH --output=logs/haps_list_fast_%A_%a.out
#SBATCH --error=logs/haps_list_fast_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-5

# FAST List-Format Haplotype Estimation for All Chromosomes
# Usage: sbatch run_list_format_all_chroms_fast.slurm <dataset> <param_file> <output_dir> [debug]
# Example: sbatch run_list_format_all_chroms_fast.slurm ZINC2 helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
# Example: sbatch run_list_format_all_chroms_fast.slurm ZINC2 helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 debug

set -euo pipefail

# Get command line arguments
if [ $# -lt 3 ] || [ $# -gt 4 ]; then
    echo "Usage: sbatch run_list_format_all_chroms_fast.slurm <dataset> <param_file> <output_dir> [debug]"
    echo "Example: sbatch run_list_format_all_chroms_fast.slurm ZINC2 helpfiles/ZINC2_haplotype_parameters.R process/ZINC2"
    echo "Example: sbatch run_list_format_all_chroms_fast.slurm ZINC2 helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 debug"
    exit 1
fi

DATASET=$1
PARAM_FILE=$2
OUTPUT_DIR=$3
DEBUG_MODE=${4:-""}

# Chromosomes to process
chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
CHR=${chrs[$((SLURM_ARRAY_TASK_ID-1))]}

mkdir -p logs

echo "[INFO] $(date) starting FAST ${CHR} for dataset ${DATASET}"
echo "[INFO] Parameter file: ${PARAM_FILE}"
echo "[INFO] Output directory: ${OUTPUT_DIR}"
echo "[INFO] Using VECTORIZED pipeline (no loops!)"
module load R/4.4.2

# 1) FAST Adaptive list-format estimation (vectorized)
if [ "$DEBUG_MODE" = "debug" ]; then
    echo "[INFO] Running in DEBUG mode - verbose output enabled"
    Rscript scripts/production/haplotype_estimation_pipeline_fast.R \
      "${CHR}" adaptive 4 "${OUTPUT_DIR}" "${PARAM_FILE}" debug
else
    echo "[INFO] Running in PRODUCTION mode - minimal output"
    Rscript scripts/production/haplotype_estimation_pipeline_fast.R \
      "${CHR}" adaptive 4 "${OUTPUT_DIR}" "${PARAM_FILE}"
fi

# 2) Smooth h4 from that adaptive list-format (writes both adapt/smooth reshaped files)
Rscript scripts/production/haplotype_estimation_pipeline_fast.R \
  smooth "${CHR}" "${PARAM_FILE}" "${OUTPUT_DIR}"

echo "[INFO] $(date) done FAST ${CHR}"
