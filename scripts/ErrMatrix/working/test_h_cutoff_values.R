#!/usr/bin/env Rscript

# Test different h_cutoff values using the existing hunk data
# This will show how h_cutoff affects error matrix behavior

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# Load the production comparison data to get fixed method results
cat("=== LOADING PRODUCTION COMPARISON DATA ===\n")
comparison_data <- readRDS("testing_positions_comparison.rds")

# Extract fixed method results for position 19780000, sample Rep01_W_F
pos_data <- comparison_data[comparison_data$pos == 19780000 & comparison_data$method == "fixed", ]
samples <- pos_data$sample[[1]]
haps_list <- pos_data$Haps[[1]]
err_list <- pos_data$Err[[1]]

if("Rep01_W_F" %in% samples) {
  idx <- which(samples == "Rep01_W_F")
  fixed_haps <- haps_list[[idx]]
  fixed_err_diag_sum <- sum(diag(err_list[[idx]]))
  cat("Fixed method results for Rep01_W_F at 19780000:\n")
  cat("Haplotypes:", paste(round(fixed_haps, 6), collapse = ", "), "\n")
  cat("Error diagonal sum:", fixed_err_diag_sum, "\n\n")
} else {
  stop("Rep01_W_F not found in fixed method data")
}

# Load the existing hunk data (h_cutoff=4)
hunk_file <- "hunk_data_chr3R_19780000_Rep01_W_F_h4.rds"
if (!file.exists(hunk_file)) {
  stop("Hunk data file not found: ", hunk_file)
}

hunk_data <- readRDS(hunk_file)
cat("Loaded hunk data with", nrow(hunk_data$df3), "rows\n\n")

# Source the est_haps_var function from the debug wrapper
source("scripts/ErrMatrix/working/debug_wrapper.R", local = TRUE)

# Test different h_cutoff values
h_cutoff_values <- c(4, 6, 8, 10)
results <- data.frame()

cat("=== TESTING DIFFERENT H_CUTOFF VALUES ===\n")

for (h_cutoff in h_cutoff_values) {
  cat("Testing h_cutoff =", h_cutoff, "\n")
  
  # Run est_haps_var with this h_cutoff value
  args_with_h_cutoff <- hunk_data$args
  args_with_h_cutoff$h_cutoff <- h_cutoff
  args_with_h_cutoff$verbose <- 1
  
  result <- do.call(est_haps_var, c(list(df3 = hunk_data$df3), args_with_h_cutoff))
  
  # Extract results
  haps <- result$Haps
  err_diag_sum <- sum(diag(result$Err))
  
  # Calculate differences from fixed method
  hap_diff <- sum(abs(haps - fixed_haps))
  err_ratio <- err_diag_sum / fixed_err_diag_sum
  
  # Store results
  results <- rbind(results, data.frame(
    h_cutoff = h_cutoff,
    hap_1 = haps[1],
    hap_2 = haps[2], 
    hap_3 = haps[3],
    hap_4 = haps[4],
    hap_5 = haps[5],
    hap_6 = haps[6],
    hap_7 = haps[7],
    hap_8 = haps[8],
    err_diag_sum = err_diag_sum,
    hap_diff_from_fixed = hap_diff,
    err_ratio_vs_fixed = err_ratio
  ))
  
  cat("  Haplotypes:", paste(round(haps, 6), collapse = ", "), "\n")
  cat("  Error diagonal sum:", err_diag_sum, "\n")
  cat("  Hap diff from fixed:", hap_diff, "\n")
  cat("  Error ratio vs fixed:", err_ratio, "\n\n")
}

# Add fixed method results for comparison
results <- rbind(results, data.frame(
  h_cutoff = "fixed",
  hap_1 = fixed_haps[1],
  hap_2 = fixed_haps[2],
  hap_3 = fixed_haps[3], 
  hap_4 = fixed_haps[4],
  hap_5 = fixed_haps[5],
  hap_6 = fixed_haps[6],
  hap_7 = fixed_haps[7],
  hap_8 = fixed_haps[8],
  err_diag_sum = fixed_err_diag_sum,
  hap_diff_from_fixed = 0,
  err_ratio_vs_fixed = 1
))

# Display results table
cat("=== RESULTS TABLE ===\n")
print(results)

# Save results
write_csv(results, "h_cutoff_debug_results.csv")
cat("\nResults saved to: h_cutoff_debug_results.csv\n")