#!/usr/bin/env Rscript

# Quick SNP Imputation Test Wrapper
# Tests imputation on 1000 SNPs starting from position 10,000,000 using fixed_50kb

suppressPackageStartupMessages(library(dplyr))

cat("=== SNP IMPUTATION TEST WRAPPER ===\n")
cat("Testing: 1000 SNPs starting from 10,000,000 using fixed_50kb estimator\n")
cat("This should run very fast (seconds, not hours)!\n\n")

# Parameters
chr <- "chr2R"
param_file <- "helpfiles/JUICE_haplotype_parameters.R" 
output_dir <- "process/JUICE/haplotype_results"
estimator <- "fixed_50kb"
start_pos <- 10000000
max_snps <- 1000

# Expected output file (testing mode creates separate file)
output_file <- file.path(output_dir, paste0("snp_imputation_", estimator, "_", chr, "_TEST.RDS"))

cat("=== INPUT/OUTPUT SUMMARY ===\n")
cat("Input haplotype file: ", file.path(output_dir, paste0(estimator, "_results_", chr, ".RDS")), "\n")
cat("Input REFALT file: process/JUICE/RefAlt.", chr, ".txt\n")
cat("Expected output file: ", output_file, "\n")
cat("Output format: RDS data frame with columns:\n")
cat("  - chr, pos, sample, observed, imputed, estimator\n")
cat("  - Will contain ~6000 rows (1000 SNPs × 6 samples)\n")
cat("  - Will show correlation between observed and imputed frequencies\n\n")

# Call the production script with testing parameters
cat("=== RUNNING COMMAND ===\n")
cmd <- paste("Rscript scripts/euchromatic_SNP_imputation_single.R", chr, param_file, output_dir, estimator, start_pos, max_snps)
cat("Command: ", cmd, "\n\n")

system2("Rscript", args = c(
  "scripts/euchromatic_SNP_imputation_single.R",
  chr,
  param_file, 
  output_dir,
  estimator,
  start_pos,
  max_snps
))

cat("\n=== POST-RUN CHECK ===\n")
if (file.exists(output_file)) {
  cat("✓ Test output file created: ", output_file, "\n")
  cat("✓ Production results safely preserved (separate test file)\n\n")
  
  # Load and summarize the results
  results <- readRDS(output_file)
  cat("✓ Results loaded: ", nrow(results), " rows\n")
  cat("✓ Columns: ", paste(names(results), collapse = ", "), "\n")
  cat("✓ Samples (", length(unique(results$sample)), "): ", paste(unique(results$sample), collapse = ", "), "\n")
  cat("✓ SNP range: ", format(min(results$pos), big.mark=","), " to ", format(max(results$pos), big.mark=","), "\n")
  cat("✓ Overall correlation: ", round(cor(results$observed, results$imputed, use = "complete.obs"), 3), "\n\n")
  
  # Show first 50 SNPs for first sample
  first_sample <- unique(results$sample)[1]
  sample_data <- results %>% 
    filter(sample == first_sample) %>% 
    arrange(pos) %>% 
    head(50)
  
  cat("=== FIRST 50 SNPs FOR SAMPLE:", first_sample, "===\n")
  cat("Position        Observed  Imputed   Diff\n")
  cat("--------        --------  -------   ----\n")
  for(i in 1:nrow(sample_data)) {
    pos_str <- format(sample_data$pos[i], width=10, big.mark=",")
    obs_str <- sprintf("%.3f", sample_data$observed[i])
    imp_str <- sprintf("%.3f", sample_data$imputed[i])
    diff_str <- sprintf("%.3f", sample_data$imputed[i] - sample_data$observed[i])
    cat(pos_str, "   ", obs_str, "   ", imp_str, "   ", diff_str, "\n")
  }
  
  # Correlation by sample
  cat("\n=== CORRELATION BY SAMPLE ===\n")
  sample_cors <- results %>%
    group_by(sample) %>%
    summarise(correlation = cor(observed, imputed, use = "complete.obs"),
              n_snps = n()) %>%
    arrange(desc(correlation))
  
  for(i in 1:nrow(sample_cors)) {
    cat(sample_cors$sample[i], ": r =", round(sample_cors$correlation[i], 3), 
        "(n =", sample_cors$n_snps[i], "SNPs)\n")
  }
  
} else {
  cat("❌ Output file not found: ", output_file, "\n")
  cat("❌ Test failed - check error messages above\n")
}
