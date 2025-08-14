#!/bin/bash
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=haplo_from_table
#SBATCH --output=logs/haplo_table_%A_%a.out
#SBATCH --error=logs/haplo_table_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-9

# Usage: sbatch scripts/haplotype_testing_from_table.sh <params_tsv_wo_header> <parfile> <outdir>
# The params TSV must have no header, with columns: chromosome<TAB>method<TAB>parameter
# Example: sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE

module load R/4.4.2

if [ $# -ne 3 ]; then
  echo "Usage: $0 <params_tsv_wo_header> <parfile> <outdir>"
  exit 1
fi

PARAMS_TSV="$1"
PARFILE="$2"
OUTDIR="$3"

mkdir -p logs

if [ ! -f "$PARAMS_TSV" ]; then
  echo "❌ Params TSV not found: $PARAMS_TSV"
  exit 1
fi
if [ ! -f "$PARFILE" ]; then
  echo "❌ Parameter file not found: $PARFILE"
  exit 1
fi

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

echo "=== Haplotype Testing (from table) ==="
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Chromosome: $CHR"
echo "Method: $METHOD"
echo "Parameter: $PARAM"
echo "Param file: $PARFILE"
echo "Outdir: $OUTDIR"
echo "Time: $(date)"

# Validate inputs
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

# Ensure input exists
REFALT_FILE="$OUTDIR/RefAlt.$CHR.txt"
if [ ! -f "$REFALT_FILE" ]; then
  echo "❌ REFALT not found: $REFALT_FILE"; exit 1
fi

# Dispatch to the appropriate single-parameter script
if [ "$METHOD" = "fixed" ]; then
  echo "Running fixed window single with ${PARAM}kb..."
  Rscript scripts/REFALT2haps.FixedWindow.Single.R "$CHR" "$PARFILE" "$OUTDIR" "$PARAM"
  STATUS=$?
  OUTFILE="$OUTDIR/fixed_window_${PARAM}kb_results_${CHR}.RDS"
else
  echo "Running adaptive window single with h_cutoff=${PARAM}..."
  Rscript scripts/REFALT2haps.AdaptWindow.Single.R "$CHR" "$PARFILE" "$OUTDIR" "$PARAM"
  STATUS=$?
  OUTFILE="$OUTDIR/adaptive_window_h${PARAM}_results_${CHR}.RDS"
fi

if [ $STATUS -eq 0 ]; then
  echo "✓ Completed successfully"
  echo "Expected output: $OUTFILE"
else
  echo "✗ Job failed"
  exit 1
fi
