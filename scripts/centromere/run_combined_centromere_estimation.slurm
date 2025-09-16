#!/bin/bash
#SBATCH --job-name=centromere_combined
#SBATCH --output=centromere_combined_%A.out
#SBATCH --error=centromere_combined_%A.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G

# Combined-Chromosome Centromere Haplotype Estimation
# Combines 2L+2R → chr2 and 3L+3R → chr3
# 
# USAGE: Run from project root directory
# sbatch centromere/run_combined_centromere_estimation.slurm

echo "=== COMBINED-CHROMOSOME CENTROMERE ESTIMATION ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# Set working directory to centromere folder
cd centromere

# Define parameters
PARAM_FILE="../helpfiles/ZINC2_haplotype_parameters.R"
INPUT_DIR="../process/ZINC2"

echo "Processing: chr2 (2L+2R), chr3 (3L+3R)"
echo ""

# Run combined-chromosome estimation
Rscript hap_estimation_combined_centromere.R $PARAM_FILE $INPUT_DIR

echo ""
echo "=== JOB COMPLETED ==="
echo "End time: $(date)"
echo "Output files:"
ls -la *.RDS
