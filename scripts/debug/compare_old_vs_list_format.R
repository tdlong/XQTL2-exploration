#!/usr/bin/env Rscript

# Compare old format vs list format haplotype results
# Check that haplotype frequencies are consistent between formats

library(dplyr)
library(readr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript compare_old_vs_list_format.R <chromosome> [position]")
}

chromosome <- args[1]
position <- if (length(args) > 1) as.numeric(args[2]) else NULL

cat("Comparing old format vs list format for chromosome:", chromosome, "\n")
if (!is.null(position)) {
  cat("Focusing on position:", position, "\n")
}

# File paths
old_adaptive_file <- paste0("process/ZINC2/haplotype_results/adaptive_window_h4_results_", chromosome, ".RDS")
old_smooth_file <- paste0("process/ZINC2/haplotype_results/smooth_h4_results_", chromosome, ".RDS")
new_adaptive_file <- paste0("process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_", chromosome, ".RDS")
new_smooth_file <- paste0("process/ZINC2/haplotype_results_list_format/smooth_h4_results_", chromosome, ".RDS")

# Check if files exist
files_to_check <- c(old_adaptive_file, old_smooth_file, new_adaptive_file, new_smooth_file)
for (file in files_to_check) {
  if (!file.exists(file)) {
    cat("WARNING: File does not exist:", file, "\n")
  }
}

# Function to extract frequencies from old format
extract_old_frequencies <- function(data, pos = NULL) {
  if (!is.null(pos)) {
    data <- data %>% dplyr::filter(pos == !!pos)
  }
  
  # Old format has columns like B1_freq, B2_freq, etc.
  freq_cols <- names(data)[grepl("_freq$", names(data))]
  
  result <- data %>%
    dplyr::select(CHROM, pos, sample, all_of(freq_cols)) %>%
    dplyr::arrange(pos, sample)
  
  return(result)
}

# Function to extract frequencies from list format
extract_list_frequencies <- function(data, pos = NULL) {
  if (!is.null(pos)) {
    data <- data %>% dplyr::filter(pos == !!pos)
  }
  
  # List format has Haps column with frequency vectors
  result <- data %>%
    dplyr::select(CHROM, pos, sample, Haps) %>%
    dplyr::mutate(
      B1_freq = purrr::map_dbl(Haps, ~ .x[1]),
      B2_freq = purrr::map_dbl(Haps, ~ .x[2]),
      B3_freq = purrr::map_dbl(Haps, ~ .x[3]),
      B4_freq = purrr::map_dbl(Haps, ~ .x[4]),
      B5_freq = purrr::map_dbl(Haps, ~ .x[5]),
      B6_freq = purrr::map_dbl(Haps, ~ .x[6]),
      B7_freq = purrr::map_dbl(Haps, ~ .x[7]),
      B8_freq = purrr::map_dbl(Haps, ~ .x[8])
    ) %>%
    dplyr::select(-Haps) %>%
    dplyr::arrange(pos, sample)
  
  return(result)
}

# Function to compare frequencies
compare_frequencies <- function(old_data, new_data, method_name) {
  cat("\n=== Comparing", method_name, "===\n")
  
  # Check dimensions
  cat("Old format dimensions:", nrow(old_data), "x", ncol(old_data), "\n")
  cat("New format dimensions:", nrow(new_data), "x", ncol(new_data), "\n")
  
  if (nrow(old_data) != nrow(new_data)) {
    cat("WARNING: Different number of rows!\n")
    return(FALSE)
  }
  
  # Get frequency columns
  freq_cols <- c("B1_freq", "B2_freq", "B3_freq", "B4_freq", 
                 "B5_freq", "B6_freq", "B7_freq", "B8_freq")
  
  # Compare each frequency column
  max_diff <- 0
  total_diff <- 0
  n_comparisons <- 0
  
  for (col in freq_cols) {
    if (col %in% names(old_data) && col %in% names(new_data)) {
      old_vals <- old_data[[col]]
      new_vals <- new_data[[col]]
      
      # Remove NA values for comparison
      valid_idx <- !is.na(old_vals) & !is.na(new_vals)
      if (sum(valid_idx) > 0) {
        diff <- abs(old_vals[valid_idx] - new_vals[valid_idx])
        max_diff <- max(max_diff, max(diff, na.rm = TRUE))
        total_diff <- total_diff + sum(diff, na.rm = TRUE)
        n_comparisons <- n_comparisons + sum(valid_idx)
      }
    }
  }
  
  cat("Maximum frequency difference:", sprintf("%.8f", max_diff), "\n")
  cat("Total frequency difference:", sprintf("%.8f", total_diff), "\n")
  cat("Number of comparisons:", n_comparisons, "\n")
  cat("Average frequency difference:", sprintf("%.8f", total_diff / max(1, n_comparisons)), "\n")
  
  # Check if frequencies are essentially identical (within machine precision)
  tolerance <- 1e-10
  if (max_diff < tolerance) {
    cat("✅ PASS: Frequencies are identical within tolerance\n")
    return(TRUE)
  } else {
    cat("❌ FAIL: Frequencies differ beyond tolerance\n")
    return(FALSE)
  }
}

# Load and compare adaptive_h4
cat("\n" , "="*50, "\n")
cat("LOADING ADAPTIVE_H4 DATA\n")
cat("="*50, "\n")

if (file.exists(old_adaptive_file) && file.exists(new_adaptive_file)) {
  old_adaptive <- readRDS(old_adaptive_file)
  new_adaptive <- readRDS(new_adaptive_file)
  
  cat("Old adaptive format structure:\n")
  str(old_adaptive)
  cat("\nNew adaptive format structure:\n")
  str(new_adaptive)
  
  # Extract frequencies
  old_adaptive_freq <- extract_old_frequencies(old_adaptive, position)
  new_adaptive_freq <- extract_list_frequencies(new_adaptive, position)
  
  # Compare
  adaptive_match <- compare_frequencies(old_adaptive_freq, new_adaptive_freq, "ADAPTIVE_H4")
} else {
  cat("Skipping adaptive_h4 comparison - files not found\n")
  adaptive_match <- NA
}

# Load and compare smooth_h4
cat("\n" , "="*50, "\n")
cat("LOADING SMOOTH_H4 DATA\n")
cat("="*50, "\n")

if (file.exists(old_smooth_file) && file.exists(new_smooth_file)) {
  old_smooth <- readRDS(old_smooth_file)
  new_smooth <- readRDS(new_smooth_file)
  
  cat("Old smooth format structure:\n")
  str(old_smooth)
  cat("\nNew smooth format structure:\n")
  str(new_smooth)
  
  # Extract frequencies
  old_smooth_freq <- extract_old_frequencies(old_smooth, position)
  new_smooth_freq <- extract_list_frequencies(new_smooth, position)
  
  # Compare
  smooth_match <- compare_frequencies(old_smooth_freq, new_smooth_freq, "SMOOTH_H4")
} else {
  cat("Skipping smooth_h4 comparison - files not found\n")
  smooth_match <- NA
}

# Summary
cat("\n" , "="*50, "\n")
cat("SUMMARY\n")
cat("="*50, "\n")

if (!is.na(adaptive_match)) {
  cat("Adaptive_h4 frequencies match:", ifelse(adaptive_match, "✅ YES", "❌ NO"), "\n")
}

if (!is.na(smooth_match)) {
  cat("Smooth_h4 frequencies match:", ifelse(smooth_match, "✅ YES", "❌ NO"), "\n")
}

if (!is.na(adaptive_match) && !is.na(smooth_match)) {
  if (adaptive_match && smooth_match) {
    cat("\n🎉 SUCCESS: All frequencies are consistent between old and new formats!\n")
  } else {
    cat("\n⚠️  WARNING: Some frequencies differ between formats!\n")
  }
}

cat("\nComparison complete.\n")
