#!/bin/bash
#SBATCH --job-name=haplotype_estimation
#SBATCH -A tdlong_lab
#SBATCH -p highmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --output=haplotype_estimation_%j.out
#SBATCH --error=haplotype_estimation_%j.err

# Haplotype Estimation Wrapper
# Runs both fixed and adaptive window haplotype estimation for a given chromosome
# 
# Usage: sbatch haplotype_estimation_wrapper.sh <chromosome> <parameter_file> <output_directory>
# Example: sbatch haplotype_estimation_wrapper.sh chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE

# Load required modules
module load R/4.4.2

# Parse command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <chromosome> <parameter_file> <output_directory>"
    echo "Example: $0 chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE"
    exit 1
fi

CHR="$1"
PARFILE="$2"
MYDIR="$3"

echo "=== Haplotype Estimation Wrapper ==="
echo "Date: $(date)"
echo "Chromosome: $CHR"
echo "Parameter file: $PARFILE"
echo "Output directory: $MYDIR"
echo ""

# Validate chromosome
valid_chromosomes=("chr2L" "chr2R" "chr3L" "chr3R" "chrX")
if [[ ! " ${valid_chromosomes[@]} " =~ " ${CHR} " ]]; then
    echo "❌ Invalid chromosome: $CHR"
    echo "Valid chromosomes: ${valid_chromosomes[*]}"
    exit 1
fi

# Check if input files exist
if [ ! -f "$MYDIR/RefAlt.$CHR.txt" ]; then
    echo "❌ Input file $MYDIR/RefAlt.$CHR.txt not found!"
    exit 1
fi

if [ ! -f "$PARFILE" ]; then
    echo "❌ Parameter file $PARFILE not found!"
    exit 1
fi

echo "✓ Input files verified"
echo ""

# Define euchromatin boundaries for each chromosome
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

echo "Euchromatin boundaries for $CHR:"
echo "  Start: ${euchromatin_start[$CHR]}"
echo "  End: ${euchromatin_end[$CHR]}"
echo ""

# Check if haplotype results already exist
FIXED_RESULTS="$MYDIR/fixed_window_results_$CHR.RDS"
ADAPTIVE_RESULTS="$MYDIR/adaptive_window_results_$CHR.RDS"

# Run Fixed Window Assessment (if results don't exist)
if [ -f "$FIXED_RESULTS" ]; then
    echo "=== Fixed Window Assessment Results Already Exist ==="
    echo "Skipping fixed window estimation: $FIXED_RESULTS"
    echo "Delete this file to force re-run"
else
    echo "=== Running Fixed Window Assessment ==="
    echo "Command: Rscript scripts/REFALT2haps.FixedWindow.R $CHR $PARFILE $MYDIR"
    Rscript scripts/REFALT2haps.FixedWindow.R $CHR $PARFILE $MYDIR
    if [ $? -eq 0 ]; then
        echo "✓ Fixed Window Assessment completed successfully"
    else
        echo "✗ Fixed Window Assessment failed"
        exit 1
    fi
fi
echo ""

# Run Adaptive Window Analysis (if results don't exist)
if [ -f "$ADAPTIVE_RESULTS" ]; then
    echo "=== Adaptive Window Analysis Results Already Exist ==="
    echo "Skipping adaptive window estimation: $ADAPTIVE_RESULTS"
    echo "Delete this file to force re-run"
else
    echo "=== Running Adaptive Window Analysis ==="
    echo "Command: Rscript scripts/REFALT2haps.AdaptWindow.R $CHR $PARFILE $MYDIR"
    Rscript scripts/REFALT2haps.AdaptWindow.R $CHR $PARFILE $MYDIR
    if [ $? -eq 0 ]; then
        echo "✓ Adaptive Window Analysis completed successfully"
    else
        echo "✗ Adaptive Window Analysis failed"
        exit 1
    fi
fi
echo ""

# Check output files
echo "=== Checking Output Files ==="
if [ -f "$MYDIR/fixed_window_results_$CHR.RDS" ]; then
    echo "✓ Fixed window results: $MYDIR/fixed_window_results_$CHR.RDS"
    echo "  Size: $(du -h $MYDIR/fixed_window_results_$CHR.RDS | cut -f1)"
else
    echo "✗ Fixed window results not found"
fi

if [ -f "$MYDIR/adaptive_window_results_$CHR.RDS" ]; then
    echo "✓ Adaptive window results: $MYDIR/adaptive_window_results_$CHR.RDS"
    echo "  Size: $(du -h $MYDIR/adaptive_window_results_$CHR.RDS | cut -f1)"
else
    echo "✗ Adaptive window results not found"
fi
echo ""

echo "=== Haplotype Estimation Complete ==="
echo "Date: $(date)"
echo "Results saved to: $MYDIR/"
echo ""
echo "Next step: Run SNP imputation with:"
echo "sbatch scripts/snp_imputation_parallel.sh"
echo ""
