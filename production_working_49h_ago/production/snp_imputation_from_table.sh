#!/bin/bash
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=snp_imputation_pipeline
#SBATCH --output=logs/snp_imputation_%A_%a.out
#SBATCH --error=logs/snp_imputation_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-10

# =============================================================================
# SNP Imputation Pipeline from Parameter Table
# =============================================================================
# 
# This script runs SNP imputation for all haplotype estimator combinations
# defined in a parameter table. It requires that haplotype estimation has
# already been completed for all estimators.
#
# WORKFLOW:
# 1. Read chromosome/method/parameter from TSV table
# 2. Verify haplotype estimation results exist
# 3. Run SNP imputation for that estimator
# 4. Output imputed SNP frequencies for euchromatic regions
#
# INPUT:
# - Parameter table (TSV): chromosome<TAB>method<TAB>parameter
# - Parameter file (R): defines founders, samples, etc.
# - Haplotype results: from haplotype estimation pipeline
# - REFALT data: allele counts from bam2bcf2REFALT.sh
#
# OUTPUT:
# - SNP imputation results: interpolated frequencies for each SNP/sample
# - Format: snp_imputation_<estimator>_<chr>.RDS
#
# USAGE:
# sbatch scripts/production/snp_imputation_from_table.sh <params_tsv> <parfile> <outdir>
#
# EXAMPLE:
# sbatch scripts/production/snp_imputation_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE
# =============================================================================

module load R/4.4.2

if [ $# -ne 3 ]; then
  echo "Usage: $0 <params_tsv> <parfile> <outdir>"
  echo ""
  echo "Arguments:"
  echo "  params_tsv - TSV file with columns: chromosome<TAB>method<TAB>parameter"
  echo "  parfile    - R parameter file defining founders, samples, etc."
  echo "  outdir     - Output directory containing haplotype_results/"
  echo ""
  echo "Example:"
  echo "  sbatch $0 helpfiles/production_slurm_params.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE"
  echo ""
  echo "Prerequisites:"
  echo "  - Haplotype estimation must be completed for all parameter combinations"
  echo "  - REFALT data files must exist in outdir/"
  exit 1
fi

PARAMS_TSV="$1"
PARFILE="$2"
OUTDIR="$3"

mkdir -p logs

# =============================================================================
# Input Validation
# =============================================================================

if [ ! -f "$PARAMS_TSV" ]; then
  echo "❌ Parameter TSV not found: $PARAMS_TSV"
  exit 1
fi
if [ ! -f "$PARFILE" ]; then
  echo "❌ Parameter file not found: $PARFILE"
  exit 1
fi
if [ ! -d "$OUTDIR" ]; then
  echo "❌ Output directory not found: $OUTDIR"
  exit 1
fi

# =============================================================================
# Parse Parameters from TSV
# =============================================================================

# Read the specific line for this array job
LINE_NUM=$SLURM_ARRAY_TASK_ID
if [ -z "$LINE_NUM" ]; then
  echo "❌ SLURM_ARRAY_TASK_ID not set - are you running via sbatch?"
  exit 1
fi

# Extract parameters from the specified line
PARAMS_LINE=$(sed -n "${LINE_NUM}p" "$PARAMS_TSV")
if [ -z "$PARAMS_LINE" ]; then
  echo "❌ No parameters found for array task $LINE_NUM"
  echo "Check that the parameter file has enough lines"
  exit 1
fi

# Parse tab-separated values
CHR=$(echo "$PARAMS_LINE" | cut -f1)
METHOD=$(echo "$PARAMS_LINE" | cut -f2)
PARAM=$(echo "$PARAMS_LINE" | cut -f3)

if [ -z "$CHR" ] || [ -z "$METHOD" ] || [ -z "$PARAM" ]; then
  echo "❌ Invalid parameter line: $PARAMS_LINE"
  echo "Expected format: chromosome<TAB>method<TAB>parameter"
  exit 1
fi

echo "=== SNP Imputation Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Parameters from: $PARAMS_TSV (line $LINE_NUM)"
echo "Chromosome: $CHR"
echo "Method: $METHOD"
echo "Parameter: $PARAM"
echo "Parameter file: $PARFILE"
echo "Output directory: $OUTDIR"
echo "Time: $(date)"
echo ""

# =============================================================================
# Set File Paths and Estimator Name
# =============================================================================

RESULTS_DIR="$OUTDIR/haplotype_results"

# Set haplotype file names and estimator (same logic as haplotype pipeline)
if [ "$METHOD" = "fixed" ]; then
  HAPLOTYPE_FILE="$RESULTS_DIR/fixed_window_${PARAM}kb_results_${CHR}.RDS"
  ESTIMATOR="fixed_${PARAM}kb"
  echo "Fixed window estimator: ${PARAM}kb"
elif [ "$METHOD" = "smooth_h4" ]; then
  HAPLOTYPE_FILE="$RESULTS_DIR/smooth_h4_results_${CHR}.RDS"
  ESTIMATOR="smooth_h4"
  echo "Smooth h4 estimator: 21-position sliding window smoothed adaptive_h4"
else
  HAPLOTYPE_FILE="$RESULTS_DIR/adaptive_window_h${PARAM}_results_${CHR}.RDS"
  ESTIMATOR="adaptive_h${PARAM}"
  echo "Adaptive window estimator: h_cutoff = ${PARAM}"
fi

echo "Expected haplotype file: $HAPLOTYPE_FILE"
echo "Estimator: $ESTIMATOR"
echo ""

# =============================================================================
# Verify Prerequisites
# =============================================================================

echo "=== Checking Prerequisites ==="

# Check if haplotype results exist
if [ ! -f "$HAPLOTYPE_FILE" ]; then
  echo "❌ Haplotype results not found: $HAPLOTYPE_FILE"
  echo "Run haplotype estimation first using haplotype_testing_from_table.sh"
  exit 1
fi
echo "✓ Haplotype results found: $HAPLOTYPE_FILE"

# Check if REFALT data exists
REFALT_FILE="$OUTDIR/RefAlt.${CHR}.txt"
if [ ! -f "$REFALT_FILE" ]; then
  echo "❌ REFALT data not found: $REFALT_FILE"
  echo "Run bam2bcf2REFALT.sh first to generate REFALT data"
  exit 1
fi
echo "✓ REFALT data found: $REFALT_FILE"

echo ""

# =============================================================================
# Run SNP Imputation
# =============================================================================

echo "=== Running SNP Imputation ==="
echo "Estimator: $ESTIMATOR"
echo "Command: Rscript scripts/production/euchromatic_SNP_imputation_single.R $CHR $PARFILE $RESULTS_DIR $ESTIMATOR"
echo ""

# Run the SNP imputation script
Rscript scripts/production/euchromatic_SNP_imputation_single.R "$CHR" "$PARFILE" "$RESULTS_DIR" "$ESTIMATOR"
STATUS=$?

# Expected output file
OUTPUT_FILE="$RESULTS_DIR/snp_imputation_${ESTIMATOR}_${CHR}.RDS"

if [ $STATUS -eq 0 ]; then
  echo "✓ SNP imputation completed successfully"
  if [ -f "$OUTPUT_FILE" ]; then
    echo "✓ Output file created: $OUTPUT_FILE"
    
    # Get file size for summary
    FILE_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
    echo "✓ File size: $FILE_SIZE"
  else
    echo "⚠️  SNP imputation completed but output file not found: $OUTPUT_FILE"
  fi
else
  echo "❌ SNP imputation failed (exit code: $STATUS)"
  exit 1
fi

echo ""

# =============================================================================
# Job Summary
# =============================================================================

echo "=== Job Complete ==="
echo "Chromosome: $CHR"
echo "Method: $METHOD"
echo "Parameter: $PARAM"
echo "Estimator: $ESTIMATOR"
echo "Output: $OUTPUT_FILE"
echo "Time: $(date)"
echo "Duration: $SECONDS seconds"
echo ""
