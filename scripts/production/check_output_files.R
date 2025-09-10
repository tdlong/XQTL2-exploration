#!/usr/bin/env Rscript

# =============================================================================
# CHECK OUTPUT FILES - DIAGNOSTIC SCRIPT
# =============================================================================
# 
# This script opens the output files from the haplotype estimation workflow
# and displays key diagnostics to verify the files are correct.
#
# USAGE:
# Rscript scripts/production/check_output_files.R <output_dir>
#
# EXAMPLES:
# Rscript scripts/production/check_output_files.R process/JUICE
# Rscript scripts/production/check_output_files.R process/ZINC2
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript check_output_files.R <output_dir>\nExample: Rscript check_output_files.R process/JUICE")
}

output_dir <- args[1]

cat("=== CHECKING OUTPUT FILES ===\n")
cat("Output directory:", output_dir, "\n")
cat("==========================================\n\n")

# Define expected file paths
list_results_dir <- file.path(output_dir, "haplotype_results_list_format")
adaptive_file <- file.path(list_results_dir, "adaptive_window_h4_results_chr2R.RDS")
smooth_original_file <- file.path(list_results_dir, "smooth_h4_results_chr2R.RDS")
smooth_reshaped_file <- file.path(list_results_dir, "smooth_h4", "R.haps.chr2R.out.rds")
adaptive_reshaped_file <- file.path(list_results_dir, "adapt_h4", "R.haps.chr2R.out.rds")

# Check which files exist
cat("Checking for output files...\n")
files_to_check <- list(
  "Adaptive (original)" = adaptive_file,
  "Smooth (original)" = smooth_original_file,
  "Smooth (reshaped)" = smooth_reshaped_file,
  "Adaptive (reshaped)" = adaptive_reshaped_file
)

for (name in names(files_to_check)) {
  file_path <- files_to_check[[name]]
  if (file.exists(file_path)) {
    file_size <- round(file.size(file_path) / 1024^2, 2)
    cat("✓", name, ":", basename(file_path), "(", file_size, "MB)\n")
  } else {
    cat("✗", name, ":", basename(file_path), "(NOT FOUND)\n")
  }
}
cat("\n")

# Load and examine the main adaptive results file
if (file.exists(adaptive_file)) {
  cat("=== ADAPTIVE RESULTS (ORIGINAL FORMAT) ===\n")
  adaptive_data <- readRDS(adaptive_file)
  
  cat("Data frame overview:\n")
  print(adaptive_data)
  cat("\n")
  
  cat("Number of rows:", nrow(adaptive_data), "\n")
  cat("Number of unique positions:", length(unique(adaptive_data$pos)), "\n")
  cat("Number of unique samples:", length(unique(adaptive_data$sample)), "\n")
  cat("Samples:", paste(unique(adaptive_data$sample), collapse = ", "), "\n\n")
  
  # Show first two positions for first sample
  first_sample <- unique(adaptive_data$sample)[1]
  first_two_positions <- head(unique(adaptive_data$pos), 2)
  
  cat("=== DETAILED ANALYSIS: First sample, first two positions ===\n")
  cat("Sample:", first_sample, "\n")
  cat("Positions:", paste(first_two_positions, collapse = ", "), "\n\n")
  
  for (pos in first_two_positions) {
    cat("--- Position:", pos, "---\n")
    
    # Get data for this position and sample
    pos_data <- adaptive_data %>%
      filter(sample == first_sample, pos == !!pos)
    
    if (nrow(pos_data) == 1) {
      # Extract the list data
      groups <- pos_data$Groups[[1]]
      haps <- pos_data$Haps[[1]]
      err <- pos_data$Err[[1]]
      names <- pos_data$Names[[1]]
      
      cat("Groups:", paste(groups, collapse = ", "), "\n")
      cat("Number of unique groups:", length(unique(groups)), "\n")
      
      cat("Haplotype frequencies:\n")
      hap_df <- data.frame(
        Founder = names,
        Frequency = round(haps, 4)
      )
      print(hap_df)
      cat("Sum of frequencies:", round(sum(haps), 6), "\n")
      
      cat("Error matrix (first 3x3):\n")
      if (is.matrix(err) && nrow(err) >= 3 && ncol(err) >= 3) {
        print(round(err[1:3, 1:3], 6))
      } else {
        cat("Error matrix too small or invalid\n")
      }
      
      cat("Error matrix dimensions:", nrow(err), "x", ncol(err), "\n")
      cat("Error matrix condition number:", round(kappa(err), 2), "\n\n")
    } else {
      cat("No data found for this position/sample combination\n\n")
    }
  }
} else {
  cat("✗ Adaptive results file not found - cannot analyze\n")
}

# Load and examine smooth results if available
if (file.exists(smooth_original_file)) {
  cat("=== SMOOTH RESULTS (ORIGINAL FORMAT) ===\n")
  smooth_data <- readRDS(smooth_original_file)
  
  cat("Data frame overview:\n")
  print(smooth_data)
  cat("\n")
  
  cat("Number of rows:", nrow(smooth_data), "\n")
  cat("Number of unique positions:", length(unique(smooth_data$pos)), "\n")
  cat("Number of unique samples:", length(unique(smooth_data$sample)), "\n\n")
  
  # Show first two positions for first sample
  first_sample <- unique(smooth_data$sample)[1]
  first_two_positions <- head(unique(smooth_data$pos), 2)
  
  cat("=== SMOOTH ANALYSIS: First sample, first two positions ===\n")
  cat("Sample:", first_sample, "\n")
  cat("Positions:", paste(first_two_positions, collapse = ", "), "\n\n")
  
  for (pos in first_two_positions) {
    cat("--- Position:", pos, "---\n")
    
    # Get data for this position and sample
    pos_data <- smooth_data %>%
      filter(sample == first_sample, pos == !!pos)
    
    if (nrow(pos_data) == 1) {
      # Extract the list data
      groups <- pos_data$Groups[[1]]
      haps <- pos_data$Haps[[1]]
      err <- pos_data$Err[[1]]
      names <- pos_data$Names[[1]]
      
      cat("Groups:", paste(groups, collapse = ", "), "\n")
      cat("Number of unique groups:", length(unique(groups)), "\n")
      
      cat("Haplotype frequencies:\n")
      hap_df <- data.frame(
        Founder = names,
        Frequency = round(haps, 4)
      )
      print(hap_df)
      cat("Sum of frequencies:", round(sum(haps), 6), "\n")
      
      cat("Error matrix (first 3x3):\n")
      if (is.matrix(err) && nrow(err) >= 3 && ncol(err) >= 3) {
        print(round(err[1:3, 1:3], 6))
      } else {
        cat("Error matrix too small or invalid\n")
      }
      
      cat("Error matrix dimensions:", nrow(err), "x", ncol(err), "\n")
      cat("Error matrix condition number:", round(kappa(err), 2), "\n\n")
    } else {
      cat("No data found for this position/sample combination\n\n")
    }
  }
} else {
  cat("✗ Smooth results file not found - cannot analyze\n")
}

# Check reshaped files if available
if (file.exists(smooth_reshaped_file)) {
  cat("=== SMOOTH RESULTS (RESHAPED FORMAT) ===\n")
  smooth_reshaped <- readRDS(smooth_reshaped_file)
  
  cat("Data frame overview:\n")
  print(smooth_reshaped)
  cat("\n")
  
  cat("Number of positions:", nrow(smooth_reshaped), "\n")
  if (nrow(smooth_reshaped) > 0) {
    cat("Samples per position:", length(smooth_reshaped$sample[[1]]), "\n")
    cat("First position samples:", paste(smooth_reshaped$sample[[1]], collapse = ", "), "\n")
  }
  cat("\n")
}

if (file.exists(adaptive_reshaped_file)) {
  cat("=== ADAPTIVE RESULTS (RESHAPED FORMAT) ===\n")
  adaptive_reshaped <- readRDS(adaptive_reshaped_file)
  
  cat("Data frame overview:\n")
  print(adaptive_reshaped)
  cat("\n")
  
  cat("Number of positions:", nrow(adaptive_reshaped), "\n")
  if (nrow(adaptive_reshaped) > 0) {
    cat("Samples per position:", length(adaptive_reshaped$sample[[1]]), "\n")
    cat("First position samples:", paste(adaptive_reshaped$sample[[1]], collapse = ", "), "\n")
  }
  cat("\n")
}

cat("=== DIAGNOSTIC COMPLETE ===\n")
cat("✓ Checked all available output files\n")
cat("✓ Displayed key diagnostics for verification\n")
