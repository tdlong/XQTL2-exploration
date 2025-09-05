#!/bin/bash

# SLURM wrapper for list-format haplotype estimation
# Runs the new list-format estimator for all positions

#SBATCH --job-name=list_haps
#SBATCH --output=logs/list_haps_%j.out
#SBATCH --error=logs/list_haps_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Usage: sbatch run_list_format_haplotype_estimation_slurm.sh <chr> <method> <parameter> <output_dir> <param_file>
# Example: sbatch run_list_format_haplotype_estimation_slurm.sh chr2R adaptive 4 process/ZINC2 helpfiles/JUICE_haplotype_parameters.R

# Check arguments
if [ $# -ne 5 ]; then
    echo "Usage: sbatch run_list_format_haplotype_estimation_slurm.sh <chr> <method> <parameter> <output_dir> <param_file>"
    echo "Example: sbatch run_list_format_haplotype_estimation_slurm.sh chr2R adaptive 4 process/ZINC2 helpfiles/JUICE_haplotype_parameters.R"
    exit 1
fi

CHR=$1
METHOD=$2
PARAMETER=$3
OUTPUT_DIR=$4
PARAM_FILE=$5

echo "=== LIST-FORMAT HAPLOTYPE ESTIMATION SLURM JOB ==="
echo "Chromosome: $CHR"
echo "Method: $METHOD"
echo "Parameter: $PARAMETER"
echo "Output directory: $OUTPUT_DIR"
echo "Parameter file: $PARAM_FILE"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# Create output directories
mkdir -p $OUTPUT_DIR/list_results
mkdir -p logs

# Load R module
module load R/4.3.0

# Set R options
export R_LIBS_USER=/home/$USER/R/x86_64-pc-linux-gnu-library/4.3

# Run haplotype estimation for all positions
echo "Running haplotype estimation for all positions on $CHR..."
echo "The R script will determine scan positions from euchromatin boundaries"
echo ""

# Run the main script in all-positions mode (R script handles boundaries internally)
Rscript scripts/production/run_haplotype_estimation_list_format.R \
    $CHR $METHOD $PARAMETER $OUTPUT_DIR $PARAM_FILE

# Check if successful
if [ $? -eq 0 ]; then
    echo "✓ All positions completed successfully"
else
    echo "✗ Haplotype estimation failed"
    exit 1
fi

echo ""
echo "=== JOB COMPLETED ==="
echo "End time: $(date)"
echo "All positions processed"
