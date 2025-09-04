#!/bin/bash
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=list_format_pipeline
#SBATCH --output=logs/list_format_pipeline_%A.out
#SBATCH --error=logs/list_format_pipeline_%A.err
#SBATCH --time=24:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2

# =============================================================================
# List Format Pipeline - SLURM Job
# =============================================================================
# 
# This script runs the complete list format pipeline for adaptive_h4 and smooth_h4:
# 1. Run adaptive_h4 haplotype estimation (resource intensive)
# 2. Generate smooth_h4 from adaptive_h4 results (efficient)
# 3. Run adaptive_h4 SNP imputation
# 4. Run smooth_h4 SNP imputation
#
# WORKFLOW:
# 1. Run adaptive_h4 haplotype estimation with list format output
# 2. Generate smooth_h4 by averaging adaptive_h4 results
# 3. Run SNP imputation for both methods sequentially
# 4. Save results to hap_list_results/ and snp_list_results/
#
# INPUT:
# - Parameter file (R): defines founders, samples, etc.
# - Observed euchromatic data: from REFALT pipeline
#
# OUTPUT:
# - Haplotype results: adaptive_h4_list_format_<chr>.RDS, smooth_h4_list_format_<chr>.RDS
# - SNP imputation results: snp_imputation_adaptive_h4_list_format_<chr>.RDS, snp_imputation_smooth_h4_list_format_<chr>.RDS
#
# USAGE:
# sbatch scripts/production/run_list_format_pipeline_slurm.sh <chr> <param_file> <output_dir>
#
# EXAMPLE:
# sbatch scripts/production/run_list_format_pipeline_slurm.sh chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
# =============================================================================

module load R/4.4.2

if [ $# -ne 3 ]; then
  echo "Usage: $0 <chr> <param_file> <output_dir>"
  echo ""
  echo "Arguments:"
  echo "  chr        - Chromosome to process (e.g., chr2R)"
  echo "  param_file - R parameter file defining founders, samples, etc."
  echo "  output_dir - Output directory (e.g., process/ZINC2)"
  echo ""
  echo "Example:"
  echo "  sbatch $0 chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2"
  exit 1
fi

CHR=$1
PARAM_FILE=$2
OUTPUT_DIR=$3

echo "=== LIST FORMAT PIPELINE - SLURM JOB ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Chromosome: $CHR"
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo ""

# Create logs directory if it doesn't exist
mkdir -p logs

# Create results directories
mkdir -p $OUTPUT_DIR/hap_list_results
mkdir -p $OUTPUT_DIR/snp_list_results

# Step 1: Run adaptive_h4 haplotype estimation (resource intensive)
echo "Step 1: Running adaptive_h4 haplotype estimation..."
echo "Start time: $(date)"
Rscript scripts/production/run_haplotype_estimation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR adaptive_h4

if [ $? -ne 0 ]; then
    echo "ERROR: adaptive_h4 haplotype estimation failed"
    echo "End time: $(date)"
    exit 1
fi

echo "✓ adaptive_h4 haplotype estimation completed"
echo "End time: $(date)"
echo ""

# Step 2: Generate smooth_h4 from adaptive_h4 results (efficient)
echo "Step 2: Generating smooth_h4 from adaptive_h4 results..."
echo "Start time: $(date)"
Rscript scripts/production/run_haplotype_estimation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR smooth_h4

if [ $? -ne 0 ]; then
    echo "ERROR: smooth_h4 generation failed"
    echo "End time: $(date)"
    exit 1
fi

echo "✓ smooth_h4 generation completed"
echo "End time: $(date)"
echo ""

# Step 3: Run adaptive_h4 SNP imputation
echo "Step 3: Running adaptive_h4 SNP imputation..."
echo "Start time: $(date)"
Rscript scripts/production/run_snp_imputation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR adaptive_h4

if [ $? -ne 0 ]; then
    echo "ERROR: adaptive_h4 SNP imputation failed"
    echo "End time: $(date)"
    exit 1
fi

echo "✓ adaptive_h4 SNP imputation completed"
echo "End time: $(date)"
echo ""

# Step 4: Run smooth_h4 SNP imputation
echo "Step 4: Running smooth_h4 SNP imputation..."
echo "Start time: $(date)"
Rscript scripts/production/run_snp_imputation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR smooth_h4

if [ $? -ne 0 ]; then
    echo "ERROR: smooth_h4 SNP imputation failed"
    echo "End time: $(date)"
    exit 1
fi

echo "✓ smooth_h4 SNP imputation completed"
echo "End time: $(date)"
echo ""

# Final summary
echo "=== LIST FORMAT PIPELINE COMPLETE ==="
echo "Job ID: $SLURM_JOB_ID"
echo "End time: $(date)"
echo ""
echo "Results saved to:"
echo "  Haplotype results: $OUTPUT_DIR/hap_list_results/"
echo "    - adaptive_h4_list_format_$CHR.RDS"
echo "    - smooth_h4_list_format_$CHR.RDS"
echo "  SNP imputation results: $OUTPUT_DIR/snp_list_results/"
echo "    - snp_imputation_adaptive_h4_list_format_$CHR.RDS"
echo "    - snp_imputation_smooth_h4_list_format_$CHR.RDS"
echo ""

# Check file sizes
echo "File sizes:"
ls -lh $OUTPUT_DIR/hap_list_results/*$CHR.RDS 2>/dev/null || echo "  No haplotype files found"
ls -lh $OUTPUT_DIR/snp_list_results/*$CHR.RDS 2>/dev/null || echo "  No SNP imputation files found"
