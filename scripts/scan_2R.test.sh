#!/bin/bash
#SBATCH --job-name=chr2R_scan
#SBATCH -A tdlong_lab
#SBATCH -p highmem
#SBATCH --cpus-per-task=3
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=scan_2R_%j.out
#SBATCH --error=scan_2R_%j.err

# Load required modules
module load R/4.4.2

# Set working directory and parameters
CHR="chr2R"
PARFILE="helpfiles/JUICE/JUICE_haplotype_parameters.R"
MYDIR="process/JUICE"

echo "=== Starting chr2R Chromosome Scan ==="
echo "Date: $(date)"
echo "Chromosome: $CHR"
echo "Parameter file: $PARFILE"
echo "Output directory: $MYDIR"
echo ""

# Check if input files exist
if [ ! -f "$MYDIR/RefAlt.$CHR.txt" ]; then
    echo "ERROR: Input file $MYDIR/RefAlt.$CHR.txt not found!"
    exit 1
fi

if [ ! -f "$PARFILE" ]; then
    echo "ERROR: Parameter file $PARFILE not found!"
    exit 1
fi

echo "✓ Input files verified"
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
    echo "Command: Rscript scripts/REFALT2haps.Assessment.R $CHR $PARFILE $MYDIR"
    Rscript scripts/REFALT2haps.Assessment.R $CHR $PARFILE $MYDIR
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

# Run Evaluation Script
echo "=== Running Haplotype Estimator Evaluation ==="
echo "Command: Rscript scripts/evaluate_haplotype_estimators.R $CHR $PARFILE $MYDIR"
Rscript scripts/evaluate_haplotype_estimators.R $CHR $PARFILE $MYDIR
if [ $? -eq 0 ]; then
    echo "✓ Evaluation completed successfully"
else
    echo "✗ Evaluation failed"
    exit 1
fi
echo ""

# Check output files
echo "=== Checking Output Files ==="
if [ -f "$MYDIR/fixed_window_results_$CHR.RDS" ]; then
    echo "✓ Fixed window results: $MYDIR/fixed_window_results_$CHR.RDS"
else
    echo "✗ Fixed window results not found"
fi

if [ -f "$MYDIR/adaptive_window_results_$CHR.RDS" ]; then
    echo "✓ Adaptive window results: $MYDIR/adaptive_window_results_$CHR.RDS"
else
    echo "✗ Adaptive window results not found"
fi

if [ -f "$MYDIR/evaluation_table_$CHR.RDS" ]; then
    echo "✓ Evaluation table: $MYDIR/evaluation_table_$CHR.RDS"
else
    echo "✗ Evaluation table not found"
fi

if [ -f "$MYDIR/mse_results_$CHR.RDS" ]; then
    echo "✓ MSE results: $MYDIR/mse_results_$CHR.RDS"
else
    echo "✗ MSE results not found"
fi
echo ""

echo "=== chr2R Chromosome Scan Complete ==="
echo "Date: $(date)"
echo "Results saved to: $MYDIR/"
echo ""

# Optional: Show file sizes
echo "=== Output File Sizes ==="
if [ -f "$MYDIR/fixed_window_results_$CHR.RDS" ]; then
    echo "Fixed window results: $(du -h $MYDIR/fixed_window_results_$CHR.RDS | cut -f1)"
fi
if [ -f "$MYDIR/adaptive_window_results_$CHR.RDS" ]; then
    echo "Adaptive window results: $(du -h $MYDIR/adaptive_window_results_$CHR.RDS | cut -f1)"
fi
