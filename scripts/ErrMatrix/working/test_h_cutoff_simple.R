#!/usr/bin/env Rscript

# Test different h_cutoff values and compare with fixed method results
# Simple approach: run extractor + debug wrapper for each h_cutoff value

suppressPackageStartupMessages({
  library(tidyverse)
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

# Test different h_cutoff values
h_cutoff_values <- c(4, 6, 8, 10)
results <- data.frame()

cat("=== TESTING DIFFERENT H_CUTOFF VALUES ===\n")

for (h_cutoff in h_cutoff_values) {
  cat("Testing h_cutoff =", h_cutoff, "\n")
  
  # Run the extractor with this h_cutoff value
  extract_cmd <- paste0("Rscript scripts/ErrMatrix/working/extract_hunk.r chr3R 19780000 Rep01_W_F ", 
                       h_cutoff, " helpfiles/ZINC2_haplotype_parameters.R process/ZINC2")
  
  cat("Running extractor...\n")
  system(extract_cmd)
  
  # Check if the hunk data file was created
  hunk_file <- paste0("hunk_data_chr3R_19780000_Rep01_W_F_h", h_cutoff, ".rds")
  
  if (file.exists(hunk_file)) {
    # Modify the debug wrapper to use this hunk file
    debug_script <- readLines("scripts/ErrMatrix/working/debug_wrapper.R")
    debug_script[388] <- paste0('hunk_file <- "', hunk_file, '"')
    writeLines(debug_script, "temp_debug_wrapper.R")
    
    # Run the debug wrapper
    cat("Running debug wrapper...\n")
    debug_output <- system("Rscript temp_debug_wrapper.R", intern = TRUE)
    
    # Parse the output to extract results
    final_result_line <- debug_output[grep("Haps:", debug_output)]
    err_diag_line <- debug_output[grep("Error matrix diagonal sum:", debug_output)]
    
    if (length(final_result_line) > 0 && length(err_diag_line) > 0) {
      # Extract haplotypes
      haps_str <- gsub(".*Haps: ", "", final_result_line)
      haps <- as.numeric(strsplit(haps_str, ", ")[[1]])
      
      # Extract error diagonal sum
      err_diag_sum <- as.numeric(gsub(".*Error matrix diagonal sum: ", "", err_diag_line))
      
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
      
    } else {
      cat("  ERROR: Could not parse debug output\n")
    }
    
    # Clean up temp file
    if (file.exists("temp_debug_wrapper.R")) {
      file.remove("temp_debug_wrapper.R")
    }
    
  } else {
    cat("  ERROR: Hunk data file not found:", hunk_file, "\n")
  }
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
write_csv(results, "h_cutoff_comparison_results.csv")
cat("\nResults saved to: h_cutoff_comparison_results.csv\n")
