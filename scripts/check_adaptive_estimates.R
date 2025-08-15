#!/usr/bin/env Rscript

# Simple script to check if adaptive window haplotype estimates are identical
# This will load the results and compare them directly

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("=== Checking Adaptive Window Haplotype Estimates ===\n")

# Load haplotype results for different adaptive methods
adaptive_methods <- c("adaptive_h4", "adaptive_h6", "adaptive_h8", "adaptive_h10")
results_list <- list()

for (method in adaptive_methods) {
  file_path <- paste0("process/JUICE/haplotype_results/", method, "_results_chr2R.RDS")
  
  if (file.exists(file_path)) {
    cat("Loading", method, "...\n")
    results <- readRDS(file_path)
    results_list[[method]] <- results
    cat("  ✓ Loaded", nrow(results), "rows\n")
  } else {
    cat("⚠️  File not found:", file_path, "\n")
  }
}

if (length(results_list) == 0) {
  stop("No adaptive window results found!")
}

cat("\n=== Comparing Results ===\n")

# Compare the first two methods
if (length(results_list) >= 2) {
  method1 <- names(results_list)[1]
  method2 <- names(results_list)[2]
  
  cat("Comparing", method1, "vs", method2, "\n")
  
  # Get founder columns
  founder_cols <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
  
  # Compare founder frequencies
  for (founder in founder_cols) {
    if (founder %in% names(results_list[[method1]]) && founder %in% names(results_list[[method2]])) {
      values1 <- results_list[[method1]][[founder]]
      values2 <- results_list[[method2]][[founder]]
      
      # Check if they're identical
      identical_check <- identical(values1, values2)
      cat("  ", founder, ":", ifelse(identical_check, "IDENTICAL", "DIFFERENT"), "\n")
      
      if (!identical_check) {
        # Show some differences
        non_na1 <- !is.na(values1)
        non_na2 <- !is.na(values2)
        both_non_na <- non_na1 & non_na2
        
        if (sum(both_non_na) > 0) {
          diff_indices <- which(abs(values1[both_non_na] - values2[both_non_na]) > 1e-10)
          if (length(diff_indices) > 0) {
            cat("    First few differences:\n")
            for (i in head(diff_indices, 5)) {
              cat("      Index", i, ":", values1[both_non_na][i], "vs", values2[both_non_na][i], "\n")
            }
          }
        }
      }
    }
  }
  
  # Check if all founder frequencies are identical
  all_identical <- TRUE
  for (founder in founder_cols) {
    if (founder %in% names(results_list[[method1]]) && founder %in% names(results_list[[method2]])) {
      if (!identical(results_list[[method1]][[founder]], results_list[[method2]][[founder]])) {
        all_identical <- FALSE
        break
      }
    }
  }
  
  cat("\nOverall result:", ifelse(all_identical, "ALL FOUNDER FREQUENCIES ARE IDENTICAL", "FOUNDER FREQUENCIES DIFFER"), "\n")
  
  if (all_identical) {
    cat("⚠️  This suggests a bug - different h_cutoff values should produce different results!\n")
  }
}

# Show summary statistics for each method
cat("\n=== Summary Statistics ===\n")

for (method in names(results_list)) {
  cat("\n", method, ":\n")
  results <- results_list[[method]]
  
  # Count non-NA estimates for each founder
  founder_cols <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
  for (founder in founder_cols) {
    if (founder %in% names(results)) {
      non_na_count <- sum(!is.na(results[[founder]]))
      total_count <- nrow(results)
      cat("  ", founder, ":", non_na_count, "/", total_count, "estimates (", round(non_na_count/total_count*100, 1), "%)\n")
    }
  }
  
  # Show h_cutoff values if present
  if ("h_cutoff" %in% names(results)) {
    unique_cutoffs <- unique(results$h_cutoff)
    cat("  h_cutoff values:", paste(unique_cutoffs, collapse = ", "), "\n")
  }
  
  # Show n_groups if present
  if ("n_groups" %in% names(results)) {
    unique_groups <- unique(results$n_groups)
    cat("  n_groups values:", paste(unique_groups, collapse = ", "), "\n")
  }
}

cat("\n=== Check Complete ===\n")
