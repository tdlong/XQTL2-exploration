#!/usr/bin/env bash

#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --job-name=compare_adapt_reshaped
#SBATCH --output=logs/compare_adapt_reshaped_%A_%a.out
#SBATCH --error=logs/compare_adapt_reshaped_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-5

set -euo pipefail

# Parse command line arguments
if [ $# -lt 1 ]; then
    echo "Usage: sbatch run_comparison_all_chroms.slurm <output_dir> [limit]"
    echo "Example: sbatch run_comparison_all_chroms.slurm process/ZINC2"
    echo "Example: sbatch run_comparison_all_chroms.slurm process/ZINC2 100"
    exit 1
fi

OUTPUT_DIR="$1"
LIMIT="${2:-}"  # Optional limit parameter

# Chromosomes to process
chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
CHR=${chrs[$((SLURM_ARRAY_TASK_ID-1))]}

mkdir -p logs

echo "[INFO] $(date) starting comparison for ${CHR}"
echo "[INFO] Using output_dir: ${OUTPUT_DIR}"
if [ -n "${LIMIT}" ]; then
    echo "[INFO] Using limit: ${LIMIT}"
fi

module load R/4.4.2

# Run comparison script
if [ -n "${LIMIT}" ]; then
    Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R "${CHR}" "${OUTPUT_DIR}" "${LIMIT}"
else
    Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R "${CHR}" "${OUTPUT_DIR}"
fi

echo "[INFO] $(date) done ${CHR}"
