#!/bin/bash
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=haplotype_pipeline
#SBATCH --output=logs/haplotype_pipeline_%A_%a.out
#SBATCH --error=logs/haplotype_pipeline_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-9

# =============================================================================
# Haplotype Testing and SNP Imputation Pipeline
# =============================================================================
# 
# This script is the main wrapper for the XQTL haplotype analysis pipeline.
# It reads a parameter table and runs haplotype estimation for each combination,
# optionally followed by SNP imputation.
#
# WORKFLOW:
# 1. Read chromosome/method/parameter from TSV table
# 2. Run haplotype estimation (fixed or adaptive window)
# 3. Optionally run SNP imputation for euchromatic regions
# 4. Output haplotype estimates and imputed SNP frequencies
#
# INPUT:
# - Parameter table (TSV): chromosome<TAB>method<TAB>parameter
# - Parameter file (R): defines founders, samples, etc.
# - REFALT data: allele counts from bam2bcf2REFALT.sh
#
# OUTPUT:
# - Haplotype estimates: founder frequencies at each position
# - SNP imputation: interpolated frequencies for each SNP
# - All results filtered to euchromatic regions only
#
# USAGE:
# sbatch scripts/haplotype_testing_from_table.sh <params_tsv> <parfile> <outdir> [run_imputation]
#
# EXAMPLES:
# # Haplotype estimation only
# sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
#
# # Haplotype estimation + SNP imputation
# sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE yes
# =============================================================================

module load R/4.4.2

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
  echo "Usage: $0 <params_tsv> <parfile> <outdir> [run_imputation]"
  echo ""
  echo "Arguments:"
  echo "  params_tsv     - TSV file with columns: chromosome<TAB>method<TAB>parameter"
  echo "  parfile        - R parameter file defining founders, samples, etc."
  echo "  outdir         - Output directory for results"
  echo "  run_imputation - Optional: 'yes' to run SNP imputation after haplotype estimation"
  echo ""
  echo "Examples:"
  echo "  # Haplotype estimation only"
  echo "  sbatch $0 helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE"
  echo ""
  echo "  # Haplotype estimation + SNP imputation"
  echo "  sbatch $0 helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE yes"
  exit 1
fi

PARAMS_TSV="$1"
PARFILE="$2"
OUTDIR="$3"
RUN_IMPUTATION="${4:-no}"

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

# =============================================================================
# Euchromatin Boundaries (Release 6 coordinates)
# =============================================================================

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

# =============================================================================
# Parse Parameters from TSV
# =============================================================================

# Map array index to row directly (no header)
LINE_NUM=$SLURM_ARRAY_TASK_ID
ROW=$(sed -n "${LINE_NUM}p" "$PARAMS_TSV")
if [ -z "$ROW" ]; then
  echo "❌ No row found for array index $SLURM_ARRAY_TASK_ID"
  exit 1
fi

CHR=$(echo "$ROW" | awk -F"\t" '{print $1}')
METHOD=$(echo "$ROW" | awk -F"\t" '{print $2}')
PARAM=$(echo "$ROW" | awk -F"\t" '{print $3}')

echo "=== Haplotype Pipeline Job ==="
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Chromosome: $CHR"
echo "Method: $METHOD"
echo "Parameter: $PARAM"
echo "Parameter file: $PARFILE"
echo "Output directory: $OUTDIR"
echo "Run SNP imputation: $RUN_IMPUTATION"
echo "Euchromatin region: ${euchromatin_start[$CHR]} - ${euchromatin_end[$CHR]} bp"
echo "Time: $(date)"
echo ""

# =============================================================================
# Validate Parameters
# =============================================================================

case "$CHR" in
  chr2L|chr2R|chr3L|chr3R|chrX) : ;;
  *) echo "❌ Invalid chromosome: $CHR"; exit 1;;
esac
case "$METHOD" in
  fixed|adaptive) : ;;
  *) echo "❌ Invalid method: $METHOD"; exit 1;;
esac
if ! [[ "$PARAM" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  echo "❌ Parameter must be numeric: $PARAM"; exit 1
fi

# =============================================================================
# Check Input Files
# =============================================================================

REFALT_FILE="$OUTDIR/RefAlt.$CHR.txt"
if [ ! -f "$REFALT_FILE" ]; then
  echo "❌ REFALT data not found: $REFALT_FILE"
  echo "Make sure you've run bam2bcf2REFALT.sh first"
  exit 1
fi

echo "✓ Input files verified"
echo ""

# =============================================================================
# Run Haplotype Estimation
# =============================================================================

echo "=== Step 1: Haplotype Estimation ==="
echo "Method: $METHOD"
echo "Parameter: $PARAM"
echo ""

if [ "$METHOD" = "fixed" ]; then
  echo "Running fixed window estimation with ${PARAM}kb window..."
  Rscript scripts/REFALT2haps.FixedWindow.Single.R "$CHR" "$PARFILE" "$OUTDIR" "$PARAM"
  STATUS=$?
  HAPLOTYPE_FILE="$OUTDIR/fixed_window_${PARAM}kb_results_${CHR}.RDS"
  ESTIMATOR="fixed_${PARAM}kb"
else
  echo "Running adaptive window estimation with h_cutoff = $PARAM..."
  Rscript scripts/REFALT2haps.AdaptWindow.Single.R "$CHR" "$PARFILE" "$OUTDIR" "$PARAM"
  STATUS=$?
  HAPLOTYPE_FILE="$OUTDIR/adaptive_window_h${PARAM}_results_${CHR}.RDS"
  ESTIMATOR="adaptive_h${PARAM}"
fi

if [ $STATUS -eq 0 ]; then
  echo "✓ Haplotype estimation completed successfully"
  echo "Output: $HAPLOTYPE_FILE"
else
  echo "✗ Haplotype estimation failed"
  exit 1
fi

echo ""

# =============================================================================
# Optional: Run SNP Imputation
# =============================================================================

if [ "$RUN_IMPUTATION" = "yes" ]; then
  echo "=== Step 2: SNP Imputation ==="
  echo "Estimator: $ESTIMATOR"
  echo "Euchromatin region: ${euchromatin_start[$CHR]} - ${euchromatin_end[$CHR]} bp"
  echo ""
  
  # Check if haplotype results exist
  if [ ! -f "$HAPLOTYPE_FILE" ]; then
    echo "❌ Haplotype results not found: $HAPLOTYPE_FILE"
    echo "Cannot run SNP imputation without haplotype estimates"
    exit 1
  fi
  
  # SNP imputation now works directly with raw REFALT files
  # No additional processed data required
  
  echo "Running SNP imputation for $ESTIMATOR..."
  Rscript scripts/euchromatic_SNP_imputation_single.R "$CHR" "$PARFILE" "$OUTDIR" "$ESTIMATOR"
  IMPUTATION_STATUS=$?
  
  if [ $IMPUTATION_STATUS -eq 0 ]; then
    echo "✓ SNP imputation completed successfully"
    echo "Output: $OUTDIR/snp_imputation_${ESTIMATOR}_${CHR}.RDS"
  else
    echo "✗ SNP imputation failed"
    exit 1
  fi
else
  echo "=== SNP Imputation Skipped ==="
  echo "To run SNP imputation, add 'yes' as the 4th argument"
  echo "Example: sbatch $0 $PARAMS_TSV $PARFILE $OUTDIR yes"
fi

echo ""
echo "=== Job Complete ==="
echo "Chromosome: $CHR"
echo "Method: $METHOD"
echo "Parameter: $PARAM"
echo "Estimator: $ESTIMATOR"
echo "SNP imputation: $RUN_IMPUTATION"
echo "Time: $(date)"
echo "Duration: $SECONDS seconds"
echo ""
