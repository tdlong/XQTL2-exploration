#!/bin/bash

# Run Haplotype Estimation and SNP Imputation with List Format
# Usage: ./scripts/production/run_list_format_pipeline.sh <chr> <param_file> <output_dir>

if [ $# -ne 3 ]; then
    echo "Usage: $0 <chr> <param_file> <output_dir>"
    echo "Example: $0 chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2"
    exit 1
fi

CHR=$1
PARAM_FILE=$2
OUTPUT_DIR=$3

echo "=== RUNNING LIST FORMAT PIPELINE ==="
echo "Chromosome: $CHR"
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Step 1: Run adaptive_h4 haplotype estimation
echo "Step 1: Running adaptive_h4 haplotype estimation..."
Rscript scripts/production/run_haplotype_estimation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR adaptive_h4

if [ $? -ne 0 ]; then
    echo "Error: adaptive_h4 haplotype estimation failed"
    exit 1
fi

echo ""

# Step 2: Run smooth_h4 haplotype estimation
echo "Step 2: Running smooth_h4 haplotype estimation..."
Rscript scripts/production/run_haplotype_estimation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR smooth_h4

if [ $? -ne 0 ]; then
    echo "Error: smooth_h4 haplotype estimation failed"
    exit 1
fi

echo ""

# Step 3: Run adaptive_h4 SNP imputation
echo "Step 3: Running adaptive_h4 SNP imputation..."
Rscript scripts/production/run_snp_imputation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR adaptive_h4

if [ $? -ne 0 ]; then
    echo "Error: adaptive_h4 SNP imputation failed"
    exit 1
fi

echo ""

# Step 4: Run smooth_h4 SNP imputation
echo "Step 4: Running smooth_h4 SNP imputation..."
Rscript scripts/production/run_snp_imputation_list_format.R $CHR $PARAM_FILE $OUTPUT_DIR smooth_h4

if [ $? -ne 0 ]; then
    echo "Error: smooth_h4 SNP imputation failed"
    exit 1
fi

echo ""
echo "=== LIST FORMAT PIPELINE COMPLETE ==="
echo "Results saved to:"
echo "  Haplotype results: $OUTPUT_DIR/hap_list_results/"
echo "  SNP imputation results: $OUTPUT_DIR/snp_list_results/"
