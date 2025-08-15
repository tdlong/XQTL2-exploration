#!/usr/bin/env Rscript

# Script to analyze adaptive window haplotype results
# This will examine what the results actually contain and whether they're identical

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("=== Analyzing Adaptive Window Haplotype Results ===\n")

# Check what haplotype result files exist
results_dir <- "process/JUICE/haplotype_results"
if (!dir.exists(results_dir)) {
  cat("❌ Results directory not found:", results_dir, "\n")
  quit(status = 1)
}

haplotype_files <- list.files(results_dir, pattern = "*_results_chr2R.RDS", full.names = TRUE)
cat("Found haplotype result files:\n")
for (file in haplotype_files) {
  cat("  ", basename(file), "\n")
}

# Focus on adaptive window results
adaptive_files <- haplotype_files[grepl("adaptive_window_h", haplotype_files)]
cat("\nAdaptive window files:\n")
for (file in adaptive_files) {
  cat("  ", basename(file), "\n")
}

if (length(adaptive_files) == 0) {
  cat("❌ No adaptive window results found!\n")
  quit(status = 1)
}

# Load and analyze each adaptive window result
for (file in adaptive_files) {
  cat("\n=== Analyzing", basename(file), "===\n")
  
  results <- readRDS(file)
  cat("  Total rows:", nrow(results), "\n")
  cat("  Columns:", paste(names(results), collapse = ", "), "\n")
  
  # Check for key columns
  founder_cols <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
  missing_founders <- founder_cols[!founder_cols %in% names(results)]
  if (length(missing_founders) > 0) {
    cat("  ⚠️  Missing founder columns:", paste(missing_founders, collapse = ", "), "\n")
  }
  
  # Show unique values for key columns
  if ("h_cutoff" %in% names(results)) {
    unique_cutoffs <- unique(results$h_cutoff)
    cat("  h_cutoff values:", paste(unique_cutoffs, collapse = ", "), "\n")
  }
  
  if ("n_groups" %in% names(results)) {
    unique_groups <- unique(results$n_groups)
    cat("  n_groups values:", paste(unique_groups, collapse = ", "), "\n")
  }
  
  if ("sample" %in% names(results)) {
    unique_samples <- unique(results$sample)
    cat("  Samples:", paste(unique_samples, collapse = ", "), "\n")
  }
  
  # Check position range
  if ("pos" %in% names(results)) {
    cat("  Position range:", min(results$pos, na.rm = TRUE), "-", max(results$pos, na.rm = TRUE), "\n")
    cat("  Number of unique positions:", length(unique(results$pos)), "\n")
  }
  
  # Analyze founder frequency estimates
  if (all(founder_cols %in% names(results))) {
    cat("  Founder frequency analysis:\n")
    
    for (founder in founder_cols) {
      values <- results[[founder]]
      non_na_count <- sum(!is.na(values))
      total_count <- length(values)
      cat("    ", founder, ":", non_na_count, "/", total_count, "estimates (", 
          round(non_na_count/total_count*100, 1), "%)\n")
      
      if (non_na_count > 0) {
        cat("      Range:", round(min(values, na.rm = TRUE), 6), "-", 
            round(max(values, na.rm = TRUE), 6), "\n")
        cat("      Mean:", round(mean(values, na.rm = TRUE), 6), "\n")
        cat("      SD:", round(sd(values, na.rm = TRUE), 6), "\n")
      }
    }
  }
}

# Compare adaptive window results
if (length(adaptive_files) >= 2) {
  cat("\n=== Comparing Adaptive Window Results ===\n")
  
  # Load first two adaptive files
  file1 <- adaptive_files[1]
  file2 <- adaptive_files[2]
  
  cat("Comparing", basename(file1), "vs", basename(file2), "\n")
  
  results1 <- readRDS(file1)
  results2 <- readRDS(file2)
  
  # Check if they have the same structure
  cat("  Same number of rows:", nrow(results1) == nrow(results2), "\n")
  cat("  Same columns:", identical(names(results1), names(results2)), "\n")
  
  # Compare founder frequencies
  if (all(founder_cols %in% names(results1)) && all(founder_cols %in% names(results2))) {
    cat("  Founder frequency comparison:\n")
    
    for (founder in founder_cols) {
      values1 <- results1[[founder]]
      values2 <- results2[[founder]]
      
      # Check if identical
      identical_check <- identical(values1, values2)
      cat("    ", founder, ":", ifelse(identical_check, "IDENTICAL", "DIFFERENT"), "\n")
      
      if (!identical_check) {
        # Show some differences
        non_na1 <- !is.na(values1)
        non_na2 <- !is.na(values2)
        both_non_na <- non_na1 & non_na2
        
        if (sum(both_non_na) > 0) {
          differences <- abs(values1[both_non_na] - values2[both_non_na])
          max_diff <- max(differences, na.rm = TRUE)
          mean_diff <- mean(differences, na.rm = TRUE)
          cat("      Max difference:", round(max_diff, 8), "\n")
          cat("      Mean difference:", round(mean_diff, 8), "\n")
          
          # Show first few differences
          diff_indices <- which(differences > 1e-10)
          if (length(diff_indices) > 0) {
            cat("      First few differences:\n")
            for (i in head(diff_indices, 3)) {
              cat("        Index", i, ":", round(values1[both_non_na][i], 6), "vs", 
                  round(values2[both_non_na][i], 6), "\n")
            }
          }
        }
      }
    }
  }
  
  # Check if all founder frequencies are identical
  all_identical <- TRUE
  for (founder in founder_cols) {
    if (founder %in% names(results1) && founder %in% names(results2)) {
      if (!identical(results1[[founder]], results2[[founder]])) {
        all_identical <- FALSE
        break
      }
    }
  }
  
  cat("\n  Overall result:", ifelse(all_identical, "ALL FOUNDER FREQUENCIES ARE IDENTICAL", "FOUNDER FREQUENCIES DIFFER"), "\n")
  
  if (all_identical) {
    cat("  ⚠️  This suggests a bug - different h_cutoff values should produce different results!\n")
  }
}

# Check for evidence of what algorithm was used
cat("\n=== Algorithm Analysis ===\n")

# Look for evidence of hierarchical clustering vs iterative approach
if (length(adaptive_files) > 0) {
  results <- readRDS(adaptive_files[1])
  
  # Check if n_groups varies (suggests hierarchical clustering was working)
  if ("n_groups" %in% names(results)) {
    unique_groups <- unique(results$n_groups)
    cat("  n_groups variation:", ifelse(length(unique_groups) > 1, "YES - suggests clustering worked", "NO - suggests clustering failed"), "\n")
    cat("  n_groups values:", paste(unique_groups, collapse = ", "), "\n")
  }
  
  # Check if h_cutoff values are correct
  if ("h_cutoff" %in% names(results)) {
    unique_cutoffs <- unique(results$h_cutoff)
    cat("  h_cutoff values found:", paste(unique_cutoffs, collapse = ", "), "\n")
  }
  
  # Check position coverage
  if ("pos" %in% names(results)) {
    total_positions <- length(unique(results$pos))
    cat("  Total positions processed:", total_positions, "\n")
    
    # Estimate if this looks like full chromosome or subset
    if (total_positions > 10000) {
      cat("  ✓ Looks like full chromosome coverage\n")
    } else {
      cat("  ⚠️  Looks like subset testing only\n")
    }
  }
}

cat("\n=== Analysis Complete ===\n")
