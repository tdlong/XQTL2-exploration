#!/bin/bash
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=haplotype_small
#SBATCH --output=logs/haplotype_pipeline_%A_%a.out
#SBATCH --error=logs/haplotype_pipeline_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2

# This script handles only the small window sizes (20kb, 50kb) with higher memory

PARAMS_TSV="$1"
PARFILE="$2"
OUTDIR="$3"
RUN_IMPUTATION="${4:-no}"

# Source the common functions
source scripts/haplotype_testing_common.sh

# Run the pipeline
run_pipeline

