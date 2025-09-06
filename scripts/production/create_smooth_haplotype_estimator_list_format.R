#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
})

# =============================================================================
# Create Smooth Haplotype Estimator - LIST FORMAT
# =============================================================================
# 
# This script creates a "smooth_h4" estimator in LIST FORMAT by applying 
# 21-position sliding window smoothing to adaptive_h4 LIST FORMAT results.
#
# KEY RULES FOR LIST FORMAT (in plain English):
# 
# GROUPS = Which founders are clustered together:
# - Groups = [1,2,3,4,5,6,7,8] = All 8 founders distinguishable (each in own group)
# - Groups = [1,1,1,1,1,1,1,1] = All 8 founders clustered together (same group)
#
# QUALITY CHECK (21-position sliding window):
# - Count how many of 21 positions have 8 distinguishable groups [1,2,3,4,5,6,7,8]
# - If ≥17 out of 21 positions are good → Quality OK
# - If <17 out of 21 positions are good → Quality NOT OK
#
# SMOOTH_H4 OUTPUT:
# - If Quality OK: Groups=[1,2,3,4,5,6,7,8], Haps=averaged frequencies, Err=averaged errors
# - If Quality NOT OK: Groups=[1,1,1,1,1,1,1,1], Haps=all NAs, Err=all NAs
#
# AVERAGING: Only average over the good positions (8 distinguishable groups), not all 21
#
# USAGE:
# Rscript scripts/production/create_smooth_haplotype_estimator_list_format.R <chr> <param_file> <output_dir>
#
# EXAMPLE:
# Rscript scripts/production/create_smooth_haplotype_estimator_list_format.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/production/create_smooth_haplotype_estimator_list_format.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

cat("=== CREATING SMOOTH HAPLOTYPE ESTIMATOR - LIST FORMAT ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
list_results_dir <- file.path(output_dir, "haplotype_results_list_format")

# Load the adaptive_h4 LIST FORMAT results
adaptive_file <- file.path(list_results_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))

if (!file.exists(adaptive_file)) {
  stop("Adaptive h4 LIST FORMAT results not found: ", adaptive_file)
}

cat("Loading adaptive_h4 LIST FORMAT results...\n")
adaptive_data <- readRDS(adaptive_file)

cat("✓ Adaptive h4 LIST FORMAT data loaded:", nrow(adaptive_data), "positions\n")
cat("Samples available:", paste(adaptive_data$sample[[1]], collapse = ", "), "\n")
cat("Founders:", paste(adaptive_data$Names[[1]][[1]], collapse = ", "), "\n\n")

# Get founder names from the data
founders <- adaptive_data$Names[[1]][[1]]

# Helper functions for smooth_h4
check_estimate_ok <- function(groups) {
  if (is.null(groups) || any(is.na(groups))) return(FALSE)
  length(unique(groups)) == 8 && all(sort(unique(groups)) == 1:8)
}

average_haps <- function(haps_list, founders) {
  if (length(haps_list) == 0) {
    return(set_names(rep(NA, length(founders)), founders))
  }
  avg_haps <- reduce(haps_list, `+`) / length(haps_list)
  
  # Only normalize if sum is not zero or very small
  if (sum(avg_haps) > 1e-10) {
    avg_haps <- avg_haps / sum(avg_haps)  # Normalize so frequencies sum to 1
  } else {
    # If sum is too small, set equal frequencies
    avg_haps <- rep(1/length(founders), length(founders))
  }
  
  names(avg_haps) <- founders
  return(avg_haps)
}

average_err <- function(err_list, founders) {
  if (length(err_list) == 0) {
    na_matrix <- matrix(NA, length(founders), length(founders))
    rownames(na_matrix) <- founders
    colnames(na_matrix) <- founders
    return(na_matrix)
  }
  
  # Check dimensions of first matrix
  first_err <- err_list[[1]]
  if (is.null(first_err) || !is.matrix(first_err)) {
    na_matrix <- matrix(NA, length(founders), length(founders))
    rownames(na_matrix) <- founders
    colnames(na_matrix) <- founders
    return(na_matrix)
  }
  
  # Ensure all matrices have same dimensions
  n_founders <- length(founders)
  if (nrow(first_err) != n_founders || ncol(first_err) != n_founders) {
    na_matrix <- matrix(NA, n_founders, n_founders)
    rownames(na_matrix) <- founders
    colnames(na_matrix) <- founders
    return(na_matrix)
  }
  
  avg_err <- reduce(err_list, `+`) / length(err_list)
  rownames(avg_err) <- founders
  colnames(avg_err) <- founders
  return(avg_err)
}

# Apply 21-position sliding window smoothing
cat("Applying 21-position sliding window smoothing using Groups information...\n")

# First add estimate_OK column
adaptive_data <- adaptive_data %>%
  mutate(estimate_OK = map_lgl(Groups, check_estimate_ok))

# Get unique samples and positions
unique_samples <- unique(adaptive_data$sample)
unique_positions <- sort(unique(adaptive_data$pos))
n_positions <- length(unique_positions)

# Process each sample separately
smooth_data <- map_dfr(unique_samples, function(sample_name) {
  cat("Processing sample:", sample_name, "\n")
  
  # Get data for this sample only
  sample_data <- adaptive_data %>% 
    filter(sample == sample_name) %>%
    arrange(pos)
  
  # Apply sliding window to this sample
  sample_data %>%
    mutate(
      # For each position, look at 21-position window
      window_quality = map_dbl(seq_len(n()), function(i) {
        start_idx <- max(1, i - 10)  # 10 positions before
        end_idx <- min(n(), i + 10)  # 10 positions after
        window_indices <- start_idx:end_idx
        
        # Count how many positions in window have estimate_OK = TRUE
        sum(sample_data$estimate_OK[window_indices], na.rm = TRUE)
      }),
      
      # Quality decision: OK if at least 17 out of 21 positions are good
      quality_ok = window_quality >= 17,
      
      # New Groups: depends on quality
      Groups = map2(Groups, quality_ok, function(groups, quality) {
        if (quality) {
          # Quality OK: All 8 founders distinguishable [1,2,3,4,5,6,7,8]
          return(1:8)
        } else {
          # Quality NOT OK: All founders clustered together [1,1,1,1,1,1,1,1]
          return(rep(1, length(founders)))
        }
      }),
      
      # New Haps: average over good positions in window
      Haps = map2(seq_len(n()), quality_ok, function(i, quality) {
        if (!quality) {
          # Quality NOT OK: Return all NAs
          return(set_names(rep(NA, length(founders)), founders))
        }
        
        start_idx <- max(1, i - 10)  # 10 positions before
        end_idx <- min(n(), i + 10)  # 10 positions after
        window_indices <- start_idx:end_idx
        
        # Get haps from good positions only
        valid_haps <- sample_data$Haps[window_indices][sample_data$estimate_OK[window_indices]]
        valid_haps <- valid_haps[map_lgl(valid_haps, ~ !any(is.na(.x)))]
        
        return(average_haps(valid_haps, founders))
      }),
      
      # New Err: average over good positions in window
      Err = map2(seq_len(n()), quality_ok, function(i, quality) {
        if (!quality) {
          # Quality NOT OK: Return all NAs
          na_matrix <- matrix(NA, length(founders), length(founders))
          rownames(na_matrix) <- founders
          colnames(na_matrix) <- founders
          return(na_matrix)
        }
        
        start_idx <- max(1, i - 10)  # 10 positions before
        end_idx <- min(n(), i + 10)  # 10 positions after
        window_indices <- start_idx:end_idx
        
        # Get error matrices from good positions only
        valid_errs <- sample_data$Err[window_indices][sample_data$estimate_OK[window_indices]]
        valid_errs <- valid_errs[map_lgl(valid_errs, ~ !any(is.na(.x)))]
        
        
        return(average_err(valid_errs, founders))
      }),
      
      # Names remain the same
      Names = map(seq_len(n()), ~ founders)
    ) %>%
    select(CHROM, pos, sample, Groups, Haps, Err, Names)  # Only keep required columns
})

# Save the smooth_h4 LIST FORMAT results
output_file <- file.path(list_results_dir, paste0("smooth_h4_results_", chr, ".RDS"))
saveRDS(smooth_data, output_file)

cat("✓ Smooth h4 LIST FORMAT results saved to:", output_file, "\n")

# Print summary statistics
cat("\n=== SMOOTH H4 LIST FORMAT SUMMARY ===\n")
cat("Total positions:", nrow(smooth_data), "\n")
cat("Positions with quality_ok = TRUE:", sum(smooth_data$quality_ok, na.rm = TRUE), "\n")
cat("Positions with quality_ok = FALSE:", sum(!smooth_data$quality_ok, na.rm = TRUE), "\n")

# Check that frequencies sum to 1 for valid positions
valid_positions <- smooth_data$quality_ok
if (sum(valid_positions, na.rm = TRUE) > 0) {
  freq_sums <- map_dbl(which(valid_positions), function(i) {
    haps <- smooth_data$Haps[[i]][[1]]
    if (any(is.na(haps))) return(NA)
    sum(haps)
  })
  
  cat("\nFrequency sum validation (valid positions only):\n")
  cat("Mean sum:", round(mean(freq_sums, na.rm = TRUE), 6), "\n")
  cat("Min sum:", round(min(freq_sums, na.rm = TRUE), 6), "\n")
  cat("Max sum:", round(max(freq_sums, na.rm = TRUE), 6), "\n")
}

cat("\n✓ Smooth h4 LIST FORMAT estimator creation complete!\n")
