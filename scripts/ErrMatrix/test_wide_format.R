#!/usr/bin/env Rscript

# =============================================================================
# TEST WRAPPER FOR WIDE FORMAT OPTIMIZATION
# =============================================================================
# 
# Simple wrapper to test the wide format optimization
# Usage: Rscript test_wide_format.R <chr> <output_dir> <param_file>

suppressPackageStartupMessages({
  library(tidyverse)
})

# Source the wide format optimization functions
source("scripts/ErrMatrix/wide_format_optimization.R")
source("scripts/ErrMatrix/est_haps_wide.R")

if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript test_wide_format.R <chr> <output_dir> <param_file>")
  }
  
  chr <- args[1]
  output_dir <- args[2]
  param_file <- args[3]
  
  cat("=== TESTING WIDE FORMAT OPTIMIZATION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n\n")
  
  # Test with debug mode (100 positions, 1 sample)
  cat("Running in DEBUG mode (100 positions × 1 sample)\n")
  
  # Run the optimized workflow
  adaptive_results <- run_adapt_h4_wide(
    chr = chr,
    method = "adaptive", 
    parameter = 4,
    output_dir = output_dir,
    param_file = param_file,
    debug = TRUE,  # Limit to 100 positions × 1 sample
    verbose = TRUE,
    debug_level = 2
  )
  
  cat("\n=== TEST COMPLETE ===\n")
  cat("✓ Processed", nrow(adaptive_results), "results\n")
  cat("✓ Results saved to output directory\n")
  
  # Show sample of results
  cat("\nSample results (first 5 rows):\n")
  print(head(adaptive_results, 5))
}
