#!/usr/bin/env Rscript

# Compare debug output with existing SNP imputation results
# This validates that our clean implementation produces the same results

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript compare_debug_vs_existing.R <chr> <estimator> <sample_index> <test_n_snps>")
}

chr <- args[1]
estimator <- args[2]
sample_index <- as.integer(args[3])
test_n_snps <- as.integer(args[4])

cat("=== Comparing Debug vs Existing Results ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample index:", sample_index, "\n")
cat("Test SNPs:", test_n_snps, "\n\n")

# Load parameter file
param_file <- "helpfiles/JUICE_haplotype_parameters.R"
source(param_file)
sample_name <- names_in_bam[sample_index]

cat("Sample name:", sample_name, "\n\n")

# Load existing results
existing_file <- file.path("process/JUICE/haplotype_results", 
                          paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
cat("Loading existing results from:", existing_file, "\n")
existing_results <- readRDS(existing_file)

# Check structure of existing results
cat("Existing results structure:\n")
print(str(existing_results))
cat("\n")

# Filter to our sample
existing_subset <- existing_results %>%
  filter(sample == sample_name) %>%
  select(POS, observed, imputed) %>%
  arrange(POS)

cat("Extracted", nrow(existing_subset), "SNPs from existing results for sample", sample_name, "\n\n")

cat("Extracted", nrow(existing_subset), "SNPs from existing results\n\n")

# Load debug results
debug_file <- file.path("process/JUICE/haplotype_results", 
                       paste0("snp_imputation_", estimator, "_", chr, "_DEBUG.RDS"))

if (file.exists(debug_file)) {
  clean_results <- readRDS(debug_file)
  
  cat("Loaded clean results from:", debug_file, "\n")
  cat("Clean results structure:\n")
  print(str(clean_results))
  cat("\n")
  
  # Get clean subset
  clean_subset <- clean_results %>%
    select(POS, observed, imputed) %>%
    arrange(POS)
  
  cat("Extracted", nrow(clean_subset), "SNPs from clean results\n\n")
  
  # Simple left_join comparison
  cat("=== DIRECT COMPARISON ===\n")
  comparison <- clean_subset %>%
    left_join(existing_subset, by = "POS", suffix = c("_clean", "_existing"))
  
  cat("Total SNPs in clean results:", nrow(clean_subset), "\n")
  cat("Total SNPs in existing results:", nrow(existing_subset), "\n")
  cat("SNPs in both:", nrow(comparison), "\n\n")
  
  # Show first few rows
  cat("First 10 comparisons:\n")
  print(head(comparison, 10))
  
  # Calculate differences
  comparison <- comparison %>%
    mutate(
      observed_diff = abs(observed_clean - observed_existing),
      imputed_diff = abs(imputed_clean - imputed_existing)
    )
  
  cat("=== COMPARISON SUMMARY ===\n")
  cat("Total SNPs compared:", nrow(comparison), "\n")
  cat("Mean observed difference:", mean(comparison$observed_diff, na.rm = TRUE), "\n")
  cat("Mean imputed difference:", mean(comparison$imputed_diff, na.rm = TRUE), "\n")
  cat("Max observed difference:", max(comparison$observed_diff, na.rm = TRUE), "\n")
  cat("Max imputed difference:", max(comparison$imputed_diff, na.rm = TRUE), "\n")
  
  # Show first few comparisons
  cat("\nFirst 10 comparisons:\n")
  print(comparison %>% select(POS, observed_clean, observed_existing, observed_diff, 
                             imputed_clean, imputed_existing, imputed_diff) %>% head(10))
  
  # Show any large differences
  large_diffs <- comparison %>%
    filter(imputed_diff > 0.01) %>%
    arrange(desc(imputed_diff))
  
  if (nrow(large_diffs) > 0) {
    cat("\nSNPs with large differences (>1%):\n")
    print(large_diffs %>% select(POS, imputed_clean, imputed_existing, imputed_diff))
  } else {
    cat("\nAll differences are small (≤1%)\n")
  }
  
  # Check if results are essentially identical
  if (all(comparison$imputed_diff < 1e-10)) {
    cat("\n✅ PERFECT MATCH: Clean implementation produces identical results!\n")
  } else if (all(comparison$imputed_diff < 0.01)) {
    cat("\n✅ GOOD MATCH: Clean implementation produces very similar results (differences < 1%)\n")
  } else {
    cat("\n⚠️  DIFFERENCES FOUND: Clean implementation produces different results\n")
  }
  
} else {
  cat("❌ Debug file not found:", debug_file, "\n")
  cat("Please run the clean implementation first:\n")
  cat("Rscript scripts/snp_imputation_clean.R", chr, "helpfiles/JUICE_haplotype_parameters.R process/JUICE", estimator, sample_index, test_n_snps, "\n")
}

cat("\n=== COMPARISON COMPLETE ===\n")
