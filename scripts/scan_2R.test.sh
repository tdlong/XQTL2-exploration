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
module load R/4.2.0

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

# Run Fixed Window Assessment (quiet mode)
echo "=== Running Fixed Window Assessment ==="
echo "Command: Rscript scripts/REFALT2haps.Assessment.R $CHR $PARFILE $MYDIR"
Rscript scripts/REFALT2haps.Assessment.R $CHR $PARFILE $MYDIR
if [ $? -eq 0 ]; then
    echo "✓ Fixed Window Assessment completed successfully"
else
    echo "✗ Fixed Window Assessment failed"
    exit 1
fi
echo ""

# Run Adaptive Window Analysis (quiet mode)
echo "=== Running Adaptive Window Analysis ==="
echo "Command: Rscript scripts/REFALT2haps.AdaptWindow.R $CHR $PARFILE $MYDIR"
Rscript scripts/REFALT2haps.AdaptWindow.R $CHR $PARFILE $MYDIR
if [ $? -eq 0 ]; then
    echo "✓ Adaptive Window Analysis completed successfully"
else
    echo "✗ Adaptive Window Analysis failed"
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
