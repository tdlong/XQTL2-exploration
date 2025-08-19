#!/usr/bin/env Rscript

# Quick SNP Imputation Test Wrapper
# Tests imputation on 1000 SNPs starting from position 10,000,000 using fixed_50kb

cat("=== SNP IMPUTATION TEST WRAPPER ===\n")
cat("Testing: 1000 SNPs starting from 10,000,000 using fixed_50kb estimator\n")
cat("This should run very fast!\n\n")

# Parameters
chr <- "chr2R"
param_file <- "helpfiles/JUICE_haplotype_parameters.R" 
output_dir <- "process/JUICE/haplotype_results"
estimator <- "fixed_50kb"
start_pos <- 10000000
max_snps <- 1000

# Call the production script with testing parameters
system2("Rscript", args = c(
  "scripts/euchromatic_SNP_imputation_single.R",
  chr,
  param_file, 
  output_dir,
  estimator,
  start_pos,
  max_snps
))

cat("\n=== TEST COMPLETE ===\n")
cat("If successful, output file: ", file.path(output_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS")), "\n")
