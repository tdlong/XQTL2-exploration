#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

# =============================================================================
# Run SNP Imputation for Smooth H4 Estimator
# =============================================================================
# 
# This script runs SNP imputation using the smooth_h4 haplotype estimates.
#
# USAGE:
# Rscript scripts/production/run_smooth_h4_imputation.R <chr> <param_file> <output_dir>
#
# EXAMPLE:
# Rscript scripts/production/run_smooth_h4_imputation.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/production/run_smooth_h4_imputation.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

cat("=== RUNNING SMOOTH H4 SNP IMPUTATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

# Check if smooth_h4 results exist
smooth_file <- file.path(results_dir, paste0("smooth_h4_results_", chr, ".RDS"))

if (!file.exists(smooth_file)) {
  stop("Smooth h4 results not found: ", smooth_file, "\n",
       "Please run create_smooth_haplotype_estimator.R first")
}

cat("✓ Smooth h4 results found\n")

# Run SNP imputation using the existing script
cat("Running SNP imputation for smooth_h4...\n")
system(paste("Rscript scripts/production/euchromatic_SNP_imputation_single.R", 
             chr, param_file, results_dir, "smooth_h4"))

cat("✓ SNP imputation complete for smooth_h4\n")
cat("Output: snp_imputation_smooth_h4_", chr, ".RDS\n", sep = "")
