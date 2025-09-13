#!/usr/bin/env bash

#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --job-name=base_var_wide_all
#SBATCH --output=logs/base_var_wide_all_%A_%a.out
#SBATCH --error=logs/base_var_wide_all_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-5

set -euo pipefail

# Parse command line arguments
if [ $# -lt 2 ]; then
    echo "Usage: sbatch run_all_chroms.slurm <param_file> <output_dir> [parameter]"
    echo "Example: sbatch run_all_chroms.slurm helpfiles/JUICE_haplotype_parameters.R process/JUICE"
    echo "Example: sbatch run_all_chroms.slurm helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 6"
    exit 1
fi

PARAM_FILE="$1"
OUTPUT_DIR="$2"
PARAMETER="${3:-4}"  # Default to 4 if not provided

# Chromosomes to process
chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
CHR=${chrs[$((SLURM_ARRAY_TASK_ID-1))]}

mkdir -p logs

echo "[INFO] $(date) starting ${CHR} with BASE_VAR_WIDE.R"
echo "[INFO] Using param_file: ${PARAM_FILE}"
echo "[INFO] Using output_dir: ${OUTPUT_DIR}"
echo "[INFO] Using parameter: ${PARAMETER}"
module load R/4.4.2

# Run BASE_VAR_WIDE.R in production mode (non-debug, non-verbose)
# This will run the complete workflow: adaptive estimation + smoothing
Rscript scripts/ErrMatrix/BASE_VAR_WIDE.R \
  "${CHR}" adaptive "${PARAMETER}" "${OUTPUT_DIR}" "${PARAM_FILE}" --nonverbose

echo "[INFO] $(date) done ${CHR}"
