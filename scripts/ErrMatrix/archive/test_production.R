#!/usr/bin/env Rscript

# =============================================================================
# TEST PRODUCTION VERSION
# =============================================================================
# 
# Simple wrapper to test the production version and time it
# Usage: Rscript test_production.R <chr> <output_dir> <param_file>

suppressPackageStartupMessages({
  library(tidyverse)
})

# Source the production functions
source("scripts/production/complete_haplotype_workflow.R")

if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript test_production.R <chr> <output_dir> <param_file>")
  }
  
  chr <- args[1]
  output_dir <- args[2]
  param_file <- args[3]
  
  cat("=== TESTING PRODUCTION VERSION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n\n")
  
  # Test with debug mode (100 positions, 1 sample)
  cat("Running in DEBUG mode (100 positions × 1 sample)\n")
  
  start_time <- Sys.time()
  
  # Run the production workflow
  adaptive_results <- run_adaptive_estimation(
    chr = chr,
    method = "adaptive", 
    parameter = 4,
    output_dir = output_dir,
    param_file = param_file,
    debug = TRUE,  # Limit to 100 positions × 1 sample
    verbose = TRUE,
    debug_level = 2
  )
  
  end_time <- Sys.time()
  total_time <- as.numeric(end_time - start_time, units = "secs")
  
  cat("\n=== PRODUCTION TEST COMPLETE ===\n")
  cat("✓ Processed", nrow(adaptive_results), "results\n")
  cat("✓ Total time:", round(total_time, 2), "seconds\n")
  cat("✓ Time per operation:", round(total_time / nrow(adaptive_results), 3), "seconds\n")
  
  # Show sample of results
  cat("\nSample results (first 5 rows):\n")
  print(head(adaptive_results, 5))
}
