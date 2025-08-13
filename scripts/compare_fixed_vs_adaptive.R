#!/usr/bin/env Rscript

# Compare Fixed vs Adaptive Window Estimators
# Check if duplicate issue only affects adaptive estimators

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"

cat("=== Comparing Fixed vs Adaptive Window Estimators ===\n")
cat("Chromosome:", chr, "\n\n")

# Check if fixed window results exist
fixed_file <- file.path(output_dir, paste0("fixed_window_results_", chr, ".RDS"))
adaptive_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))

cat("=== File Check ===\n")
cat("Fixed window file exists:", file.exists(fixed_file), "\n")
cat("Adaptive window file exists:", file.exists(adaptive_file), "\n\n")

# Load adaptive window results
if (file.exists(adaptive_file)) {
  adaptive_results <- read_rds(adaptive_file)
  cat("✓ Adaptive window results loaded:", nrow(adaptive_results), "rows\n")
  
  # Check structure
  cat("Adaptive columns:", paste(names(adaptive_results), collapse = ", "), "\n")
  cat("Adaptive h_cutoffs:", paste(sort(unique(adaptive_results$h_cutoff)), collapse = ", "), "\n")
  cat("Adaptive samples:", paste(sort(unique(adaptive_results$sample)), collapse = ", "), "\n\n")
  
  # Check for duplicates in adaptive
  cat("=== Adaptive Window Duplicate Analysis ===\n")
  adaptive_duplicates <- adaptive_results %>%
    group_by(h_cutoff, sample, pos, founder) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)
  
  cat("Total duplicate positions in adaptive:", nrow(adaptive_duplicates), "\n")
  if (nrow(adaptive_duplicates) > 0) {
    cat("Max duplicates per founder:", max(adaptive_duplicates$n), "\n")
    cat("Duplicate summary by h_cutoff:\n")
    print(adaptive_duplicates %>% 
            group_by(h_cutoff) %>% 
            summarise(duplicates = n(), .groups = "drop"))
  }
} else {
  cat("❌ Adaptive window file not found\n")
  adaptive_results <- NULL
}

# Load fixed window results
if (file.exists(fixed_file)) {
  fixed_results <- read_rds(fixed_file)
  cat("\n✓ Fixed window results loaded:", nrow(fixed_results), "rows\n")
  
  # Check structure
  cat("Fixed columns:", paste(names(fixed_results), collapse = ", "), "\n")
  cat("Fixed window sizes:", paste(sort(unique(fixed_results$window_size)), collapse = ", "), "\n")
  cat("Fixed samples:", paste(sort(unique(fixed_results$sample)), collapse = ", "), "\n\n")
  
  # Check for duplicates in fixed
  cat("=== Fixed Window Duplicate Analysis ===\n")
  fixed_duplicates <- fixed_results %>%
    group_by(window_size, sample, pos, founder) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)
  
  cat("Total duplicate positions in fixed:", nrow(fixed_duplicates), "\n")
  if (nrow(fixed_duplicates) > 0) {
    cat("Max duplicates per founder:", max(fixed_duplicates$n), "\n")
    cat("Duplicate summary by window size:\n")
    print(fixed_duplicates %>% 
            group_by(window_size) %>% 
            summarise(duplicates = n(), .groups = "drop"))
  } else {
    cat("✓ No duplicates found in fixed window results\n")
  }
} else {
  cat("❌ Fixed window file not found\n")
  fixed_results <- NULL
}

# Compare sample coverage
if (!is.null(adaptive_results) && !is.null(fixed_results)) {
  cat("\n=== Sample Coverage Comparison ===\n")
  
  # Check which samples are in each dataset
  adaptive_samples <- unique(adaptive_results$sample)
  fixed_samples <- unique(fixed_results$sample)
  
  cat("Samples in adaptive:", length(adaptive_samples), "\n")
  cat("Samples in fixed:", length(fixed_samples), "\n")
  
  # Check overlap
  common_samples <- intersect(adaptive_samples, fixed_samples)
  cat("Common samples:", length(common_samples), "\n")
  
  if (length(common_samples) > 0) {
    cat("Common sample:", common_samples[1], "\n")
    
    # Compare positions for this sample
    sample_name <- common_samples[1]
    
    adaptive_sample_data <- adaptive_results %>% 
      filter(sample == sample_name, h_cutoff == 6)  # Use h6 for comparison
    
    fixed_sample_data <- fixed_results %>% 
      filter(sample == sample_name, window_size == 10000)  # Use 10kb for comparison
    
    cat("\nPosition comparison for sample", sample_name, ":\n")
    cat("Adaptive h6 positions:", length(unique(adaptive_sample_data$pos)), "\n")
    cat("Fixed 10kb positions:", length(unique(fixed_sample_data$pos)), "\n")
    
    # Check for duplicates in this sample
    adaptive_sample_dups <- adaptive_sample_data %>%
      group_by(pos, founder) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1)
    
    fixed_sample_dups <- fixed_sample_data %>%
      group_by(pos, founder) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1)
    
    cat("Duplicates in adaptive sample:", nrow(adaptive_sample_dups), "\n")
    cat("Duplicates in fixed sample:", nrow(fixed_sample_dups), "\n")
  }
}

# Summary
cat("\n=== Summary ===\n")
if (!is.null(adaptive_results) && !is.null(fixed_results)) {
  if (nrow(adaptive_duplicates) > 0 && nrow(fixed_duplicates) == 0) {
    cat("✓ CONFIRMED: Duplicate issue ONLY affects adaptive window estimators\n")
    cat("  This supports the theory that it's caused by multiple window size iterations\n")
  } else if (nrow(adaptive_duplicates) > 0 && nrow(fixed_duplicates) > 0) {
    cat("⚠️  Duplicates found in BOTH estimators - this suggests a different issue\n")
  } else {
    cat("✓ No duplicates found in either estimator\n")
  }
} else {
  cat("Cannot compare - missing one or both result files\n")
}

cat("\n=== Next Steps ===\n")
cat("1. If duplicates only in adaptive: Fix the window iteration logic\n")
cat("2. If duplicates in both: Look for a different root cause\n")
cat("3. Check the haplotype calling code for the root issue\n")
