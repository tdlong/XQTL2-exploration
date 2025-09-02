#!/usr/bin/env Rscript

library(dplyr)

# Check structure of both files
cat("=== DEBUG FILE STRUCTURE ===\n")
debug_data <- readRDS("process/JUICE/haplotype_results/snp_imputation_adaptive_h8_chr2R_DEBUG.RDS")
print(str(debug_data))
cat("\nColumn names:", names(debug_data), "\n\n")

cat("=== EXISTING FILE STRUCTURE ===\n")
existing_data <- readRDS("process/JUICE/haplotype_results/snp_imputation_adaptive_h8_chr2R.RDS")
print(str(existing_data))
cat("\nColumn names:", names(existing_data), "\n\n")

# Show first few rows of each
cat("=== DEBUG FILE FIRST 3 ROWS ===\n")
print(head(debug_data, 3))

cat("\n=== EXISTING FILE FIRST 3 ROWS ===\n")
print(head(existing_data, 3))
