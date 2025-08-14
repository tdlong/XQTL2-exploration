#!/bin/bash
#SBATCH -A tdlong_lab     
#SBATCH -p highmem
#SBATCH --job-name=comprehensive_haplotype_testing
#SBATCH --output=logs/comprehensive_haplotype_%A_%a.out
#SBATCH --error=logs/comprehensive_haplotype_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-60

# =============================================================================
# Comprehensive Haplotype Testing Array Job
# =============================================================================
# 
# This script systematically tests all combinations of:
# - Chromosomes: chr2L, chr2R, chr3L, chr3R, chrX
# - Methods: fixed, adaptive  
# - Parameters: 
#   * Fixed: 10, 20, 50, 100, 200, 500 (kb)
#   * Adaptive: 2, 4, 6, 8, 10, 20 (h_cutoff)
#
# Total: 5 chromosomes × 12 estimators = 60 array jobs
#
# Usage: sbatch comprehensive_haplotype_testing.sh <parameter_file> <output_directory>
# Example: sbatch comprehensive_haplotype_testing.sh helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
#
# Output: Each job produces haplotype estimates for one specific combination
# =============================================================================

# Load modules
module load R/4.4.2

# Parse command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <parameter_file> <output_directory>"
    echo "Example: $0 helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE"
    exit 1
fi

PARFILE="$1"
MYDIR="$2"

# Create logs directory
mkdir -p logs

# =============================================================================
# Define the comprehensive testing matrix
# =============================================================================

# Chromosomes
chromosomes=("chr2L" "chr2R" "chr3L" "chr3R" "chrX")

# Euchromatin boundaries for each chromosome (release 6 coordinates)
declare -A euchromatin_start
declare -A euchromatin_end

euchromatin_start["chr2L"]=82455
euchromatin_end["chr2L"]=22011009
euchromatin_start["chr2R"]=5398184
euchromatin_end["chr2R"]=24684540
euchromatin_start["chr3L"]=158639
euchromatin_end["chr3L"]=22962476
euchromatin_start["chr3R"]=4552934
euchromatin_end["chr3R"]=31845060
euchromatin_start["chrX"]=277911
euchromatin_end["chrX"]=22628490

# Methods and parameters
fixed_params=(10 20 50 100 200 500)    # kb
adaptive_params=(2 4 6 8 10 20)        # h_cutoff

# =============================================================================
# Map array job index to specific combination
# =============================================================================

# Calculate which combination this job should run
# Array indices 1-60 map to:
# - Jobs 1-30: Fixed window (5 chromosomes × 6 parameters)
# - Jobs 31-60: Adaptive window (5 chromosomes × 6 parameters)

job_index=$SLURM_ARRAY_TASK_ID

if [ $job_index -le 30 ]; then
    # Fixed window jobs (1-30)
    method="fixed"
    # Calculate chromosome and parameter indices
    # Each chromosome gets 6 parameters, so:
    # chr_index = (job_index - 1) / 6
    # param_index = (job_index - 1) % 6
    chr_index=$(( (job_index - 1) / 6 ))
    param_index=$(( (job_index - 1) % 6 ))
    chr=${chromosomes[$chr_index]}
    param=${fixed_params[$param_index]}
    estimator="fixed_${param}kb"
else
    # Adaptive window jobs (31-60)
    method="adaptive"
    # Calculate chromosome and parameter indices
    # Each chromosome gets 6 parameters, so:
    # chr_index = (job_index - 31) / 6
    # param_index = (job_index - 31) % 6
    chr_index=$(( (job_index - 31) / 6 ))
    param_index=$(( (job_index - 31) % 6 ))
    chr=${chromosomes[$chr_index]}
    param=${adaptive_params[$param_index]}
    estimator="adaptive_h${param}"
fi

# =============================================================================
# Job Information
# =============================================================================

echo "=== Comprehensive Haplotype Testing Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Chromosome: $chr"
echo "Method: $method"
echo "Parameter: $param"
echo "Estimator: $estimator"
echo "Parameter file: $PARFILE"
echo "Output directory: $MYDIR"
echo "Time: $(date)"
echo ""

# =============================================================================
# Input Validation
# =============================================================================

# Check parameter file
if [ ! -f "$PARFILE" ]; then
    echo "❌ Parameter file not found: $PARFILE"
    exit 1
fi

# Check REFALT data
refalt_file="$MYDIR/RefAlt.$chr.txt"
if [ ! -f "$refalt_file" ]; then
    echo "❌ REFALT data not found: $refalt_file"
    exit 1
fi

echo "✓ Input files verified"
echo "Euchromatin boundaries for $chr: ${euchromatin_start[$chr]} - ${euchromatin_end[$chr]}"
echo ""

# =============================================================================
# Run Haplotype Estimation
# =============================================================================

echo "=== Running Haplotype Estimation ==="
echo "Method: $method"
echo "Parameter: $param"
echo ""

if [ "$method" = "fixed" ]; then
    # Fixed window estimation
    echo "Running fixed window estimation with ${param}kb window..."
    Rscript scripts/REFALT2haps.FixedWindow.Single.R $chr $PARFILE $MYDIR $param
    
    if [ $? -eq 0 ]; then
        echo "✓ Fixed window estimation completed successfully"
        echo "Output: $MYDIR/fixed_window_${param}kb_results_${chr}.RDS"
    else
        echo "✗ Fixed window estimation failed"
        exit 1
    fi
else
    # Adaptive window estimation
    echo "Running adaptive window estimation with h_cutoff = $param..."
    Rscript scripts/REFALT2haps.AdaptWindow.Single.R $chr $PARFILE $MYDIR $param
    
    if [ $? -eq 0 ]; then
        echo "✓ Adaptive window estimation completed successfully"
        echo "Output: $MYDIR/adaptive_window_h${param}_results_${chr}.RDS"
    else
        echo "✗ Adaptive window estimation failed"
        exit 1
    fi
fi

echo ""
echo "=== Job Complete ==="
echo "Chromosome: $chr"
echo "Method: $method"
echo "Parameter: $param"
echo "Estimator: $estimator"
echo "Time: $(date)"
echo "Duration: $SECONDS seconds"
echo ""

# =============================================================================
# Optional: Run SNP Imputation for this specific estimator
# =============================================================================
# Uncomment the following section if you want to run SNP imputation immediately
# after haplotype estimation for each combination

# echo "=== Running SNP Imputation for $estimator ==="
# Rscript scripts/euchromatic_SNP_imputation_single.R $chr $PARFILE $MYDIR $estimator
# if [ $? -eq 0 ]; then
#     echo "✓ SNP Imputation completed successfully"
#     echo "Output: $MYDIR/snp_imputation_${estimator}_${chr}.RDS"
# else
#     echo "✗ SNP Imputation failed"
# fi
