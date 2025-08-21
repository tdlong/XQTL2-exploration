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

# Filter to our sample and get the same SNPs as debug run
set.seed(1201)  # Same seed as debug mode
euchromatic_start <- 5398184
euchromatic_end <- 24684540

# Get valid SNP positions from existing results
valid_positions <- existing_results %>%
  filter(sample == sample_name, 
         POS >= euchromatic_start, 
         POS <= euchromatic_end) %>%
  pull(POS) %>%
  unique()

# Sample the same positions as debug mode
if (length(valid_positions) >= test_n_snps) {
  test_positions <- sample(valid_positions, test_n_snps)
} else {
  test_positions <- valid_positions
  cat("Warning: Only", length(valid_positions), "positions available, using all\n")
}

# Extract existing results for comparison
existing_subset <- existing_results %>%
  filter(sample == sample_name, POS %in% test_positions) %>%
  select(POS, observed, imputed) %>%
  arrange(POS)

cat("Extracted", nrow(existing_subset), "SNPs from existing results\n\n")

# Run our clean implementation and capture output
cat("Running clean implementation...\n")
debug_file <- file.path("process/JUICE/haplotype_results", 
                       paste0("snp_imputation_", estimator, "_", chr, "_DEBUG.RDS"))

if (file.exists(debug_file)) {
  clean_results <- readRDS(debug_file)
  
  # Filter to same positions
  clean_subset <- clean_results %>%
    filter(POS %in% test_positions) %>%
    select(POS, observed, imputed) %>%
    arrange(POS)
  
  cat("Loaded clean results from:", debug_file, "\n")
  cat("Extracted", nrow(clean_subset), "SNPs from clean results\n\n")
  
  # Compare results
  comparison <- clean_subset %>%
    left_join(existing_subset, by = "POS", suffix = c("_clean", "_existing")) %>%
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
