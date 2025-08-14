#!/bin/bash
#SBATCH -A tdlong_lab     
#SBATCH -p standard
#SBATCH --job-name=snp_imputation_parallel
#SBATCH --output=logs/snp_imputation_%A_%a.out
#SBATCH --error=logs/snp_imputation_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-12

# SNP Imputation Parallel Script
# Parallelizes over 12 haplotype estimators:
# 1-6: Fixed window sizes (10kb, 20kb, 50kb, 100kb, 200kb, 500kb)
# 7-12: Adaptive h_cutoffs (2, 4, 6, 8, 10, 20)

# Load modules
module load R/4.4.2

# Parse command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <chromosome> <parameter_file> <output_directory>"
    echo "Example: $0 chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE"
    exit 1
fi

# Set variables
CHR="$1"
PARFILE="$2"
MYDIR="$3"
SCRIPT_DIR="scripts"

# Create logs directory
mkdir -p logs

# Define the 12 haplotype estimators
declare -a estimators=(
    "fixed_10kb"
    "fixed_20kb" 
    "fixed_50kb"
    "fixed_100kb"
    "fixed_200kb"
    "fixed_500kb"
    "adaptive_h2"
    "adaptive_h4"
    "adaptive_h6"
    "adaptive_h8"
    "adaptive_h10"
    "adaptive_h20"
)

# Get the estimator for this job
ESTIMATOR=${estimators[$SLURM_ARRAY_TASK_ID-1]}

echo "=== SNP Imputation Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Estimator: $ESTIMATOR"
echo "Chromosome: $CHR"
echo "Parameter file: $PARFILE"
echo "Output directory: $MYDIR"
echo "Time: $(date)"
echo ""

# Check if haplotype results exist
if [[ $ESTIMATOR == fixed_* ]]; then
    # Fixed window estimator
    window_size=${ESTIMATOR#fixed_}
    haplotype_file="${MYDIR}/fixed_window_results_${CHR}.RDS"
    if [ ! -f "$haplotype_file" ]; then
        echo "❌ Fixed window haplotype results not found: $haplotype_file"
        exit 1
    fi
    echo "✓ Fixed window haplotype results found: $haplotype_file"
else
    # Adaptive window estimator
    h_cutoff=${ESTIMATOR#adaptive_h}
    haplotype_file="${MYDIR}/adaptive_window_results_${CHR}.RDS"
    if [ ! -f "$haplotype_file" ]; then
        echo "❌ Adaptive window haplotype results not found: $haplotype_file"
        exit 1
    fi
    echo "✓ Adaptive window haplotype results found: $haplotype_file"
fi

# Check if REFALT data exists
refalt_file="${MYDIR}/df3.${CHR}.RDS"
if [ ! -f "$refalt_file" ]; then
    echo "❌ REFALT data not found: $refalt_file"
    exit 1
fi
echo "✓ REFALT data found: $refalt_file"

# Run SNP imputation for this specific estimator
echo "=== Running SNP Imputation for $ESTIMATOR ==="
echo "Command: Rscript $SCRIPT_DIR/euchromatic_SNP_imputation_single.R $CHR $PARFILE $MYDIR $ESTIMATOR"

Rscript $SCRIPT_DIR/euchromatic_SNP_imputation_single.R $CHR $PARFILE $MYDIR $ESTIMATOR

if [ $? -eq 0 ]; then
    echo "✓ SNP Imputation for $ESTIMATOR completed successfully"
    echo "Output: ${MYDIR}/snp_imputation_${ESTIMATOR}_${CHR}.RDS"
else
    echo "✗ SNP Imputation for $ESTIMATOR failed"
    exit 1
fi

echo ""
echo "=== Job Complete ==="
echo "Estimator: $ESTIMATOR"
echo "Time: $(date)"
echo "Duration: $SECONDS seconds"
