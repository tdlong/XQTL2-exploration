#!/usr/bin/env Rscript

# =============================================================================
# BENCHMARK: WIDE FORMAT vs PRODUCTION
# =============================================================================
# 
# This script runs both the production and wide format versions
# on the same data and compares their performance
# 
# Usage: Rscript benchmark_wide_vs_production.R <chr> <output_dir> <param_file>

suppressPackageStartupMessages({
  library(tidyverse)
  library(microbenchmark)
})

# Source both versions
source("scripts/production/complete_haplotype_workflow.R")
source("scripts/ErrMatrix/wide_format_optimization.R")
source("scripts/ErrMatrix/est_haps_wide.R")

if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript benchmark_wide_vs_production.R <chr> <output_dir> <param_file>")
  }
  
  chr <- args[1]
  output_dir <- args[2]
  param_file <- args[3]
  
  cat("=== BENCHMARKING WIDE FORMAT vs PRODUCTION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n\n")
  
  # Load parameters
  source(param_file)
  
  # Load RefAlt data in wide format for comparison
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  cat("Loading data...\n")
  wide_data <- read_RefAlt_wide(refalt_file, founders)
  df3_wide <- wide_data$df3_wide
  founders <- wide_data$founders
  samples <- wide_data$samples
  
  # Convert to long format for production version
  df3_long <- df3_wide %>%
    dplyr::select(POS, dplyr::all_of(founders), dplyr::all_of(samples)) %>%
    tidyr::pivot_longer(cols = c(dplyr::all_of(founders), dplyr::all_of(samples)), 
                       names_to = "name", values_to = "freq") %>%
    dplyr::arrange(POS)
  
  # Define test positions (limited for benchmarking)
  euchromatin_boundaries <- list(
    chr2L = c(82455, 22011009),
    chr2R = c(5398184, 24684540),
    chr3L = c(158639, 22962476),
    chr3R = c(4552934, 31845060),
    chrX = c(277911, 22628490)
  )
  
  euchrom_start <- euchromatin_boundaries[[chr]][1]
  euchrom_end <- euchromatin_boundaries[[chr]][2]
  
  scan_start <- ceiling(euchrom_start / step) * step
  scan_end <- floor(euchrom_end / step) * step
  all_positions <- seq(scan_start, scan_end, by = step)
  
  # Limit to first 50 positions for benchmarking
  test_positions <- head(all_positions, 50)
  test_samples <- head(samples, 3)  # Test with 3 samples
  
  cat("Testing with", length(test_positions), "positions ×", length(test_samples), "samples\n")
  cat("Total operations:", length(test_positions) * length(test_samples), "\n\n")
  
  # Function to run production version
  run_production <- function() {
    results <- map_dfr(test_positions, function(test_pos) {
      map_dfr(test_samples, function(sample_name) {
        result <- estimate_haplotypes_list_format(
          pos = test_pos,
          sample_name = sample_name,
          df3 = df3_long,
          founders = founders,
          h_cutoff = 4,
          method = "adaptive",
          verbose = 0
        )
        
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = sample_name,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      })
    })
    return(results)
  }
  
  # Function to run wide format version
  run_wide_format <- function() {
    results <- map_dfr(test_positions, function(test_pos) {
      # Subset 500kb window for this position
      window_start <- max(1, test_pos - 250000)
      window_end <- min(max(df3_wide$POS), test_pos + 250000)
      
      window_data <- df3_wide %>%
        dplyr::filter(POS >= window_start & POS <= window_end)
      
      if (nrow(window_data) < 10) {
        # No data in window - return NA results for all samples
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = test_samples,
          Groups = map(test_samples, ~ rep(1, length(founders))),
          Haps = map(test_samples, ~ set_names(rep(NA, length(founders)), founders)),
          Err = map(test_samples, ~ matrix(NA, length(founders), length(founders))),
          Names = map(test_samples, ~ founders)
        ))
      }
      
      # Process all samples for this position
      sample_results <- map_dfr(test_samples, function(sample_name) {
        result <- est_haps_wide(
          pos = test_pos,
          sample_name = sample_name,
          df3_wide = window_data,
          founders = founders,
          h_cutoff = 4,
          verbose = 0
        )
        
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = sample_name,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      })
      
      return(sample_results)
    })
    return(results)
  }
  
  # Run benchmark
  cat("Running benchmark...\n")
  benchmark_results <- microbenchmark(
    production = run_production(),
    wide_format = run_wide_format(),
    times = 3,  # Run each version 3 times
    unit = "s"
  )
  
  # Display results
  cat("\n=== BENCHMARK RESULTS ===\n")
  print(benchmark_results)
  
  # Calculate speedup
  prod_times <- benchmark_results$time[benchmark_results$expr == "production"]
  wide_times <- benchmark_results$time[benchmark_results$expr == "wide_format"]
  
  prod_median <- median(prod_times) / 1e9  # Convert to seconds
  wide_median <- median(wide_times) / 1e9  # Convert to seconds
  
  speedup <- prod_median / wide_median
  
  cat("\n=== PERFORMANCE SUMMARY ===\n")
  cat("Production median time:", round(prod_median, 2), "seconds\n")
  cat("Wide format median time:", round(wide_median, 2), "seconds\n")
  cat("Speedup:", round(speedup, 2), "x\n")
  
  if (speedup > 1) {
    cat("✓ Wide format is", round(speedup, 2), "x FASTER\n")
  } else {
    cat("⚠ Wide format is", round(1/speedup, 2), "x SLOWER\n")
  }
  
  # Verify results are similar
  cat("\n=== VERIFYING RESULTS ===\n")
  prod_results <- run_production()
  wide_results <- run_wide_format()
  
  # Compare a few key metrics
  cat("Production results:", nrow(prod_results), "rows\n")
  cat("Wide format results:", nrow(wide_results), "rows\n")
  
  # Check if results have same structure
  if (nrow(prod_results) == nrow(wide_results)) {
    cat("✓ Same number of results\n")
  } else {
    cat("⚠ Different number of results\n")
  }
  
  cat("\n=== BENCHMARK COMPLETE ===\n")
}
