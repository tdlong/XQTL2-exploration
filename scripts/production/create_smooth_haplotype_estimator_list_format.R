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
# KEY RULES FOR LIST FORMAT:
# - Use Groups information instead of estimate_OK
# - Quality: ≥17 out of 21 positions must have 8 distinguishable groups (1:8)
# - Only average over positions with 8 groups (not all 21 positions)
# - Groups for smooth_h4: 1:8 if quality OK, else all 1s (fallback)
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

# Apply 21-position sliding window smoothing
cat("Applying 21-position sliding window smoothing using Groups information...\n")

smooth_data <- adaptive_data %>%
  arrange(pos) %>%
  mutate(
    # Calculate quality count: how many positions in 21-position window have 8 groups?
    quality_count = map_dbl(seq_len(n()), function(i) {
      start_idx <- max(1, i - 10)  # 10 positions before
      end_idx <- min(n(), i + 10)  # 10 positions after
      window_positions <- start_idx:end_idx
      
      # Check if each position has 8 distinguishable groups (1:8)
      positions_with_8_groups <- map_lgl(window_positions, function(j) {
        groups_at_pos <- adaptive_data$Groups[[j]][[1]]
        if (is.null(groups_at_pos) || any(is.na(groups_at_pos))) return(FALSE)
        # Check if all 8 founders are in different groups (1:8)
        length(unique(groups_at_pos)) == 8 && all(sort(unique(groups_at_pos)) == 1:8)
      })
      
      sum(positions_with_8_groups, na.rm = TRUE)
    }),
    
    # New quality: OK if at least 17 out of 21 positions are OK
    quality_ok = quality_count >= 17,
    
    # For smooth_h4, groups depend on the quality
    Groups = list(map(seq_along(sample[[1]]), function(i) {
      if (quality_ok) {
        # If smooth estimate is OK, use 1:8 (all founders distinguishable)
        return(1:8)
      } else {
        # If smooth estimate is not OK, use all 1s (fallback)
        return(rep(1, length(founders)))
      }
    })),
    
    # Average the haplotype frequencies using 21-position sliding window
    Haps = list(map(seq_along(sample[[1]]), function(i) {
      start_idx <- max(1, i - 10)  # 10 positions before
      end_idx <- min(nrow(adaptive_data), i + 10)  # 10 positions after
      
      # Get haplotype estimates for this sample across the 21-position window
      window_haps <- map(start_idx:end_idx, function(j) {
        adaptive_data$Haps[[j]][[1]][[i]]
      })
      
      # Only use positions where the original estimate had 8 groups
      valid_positions <- map_lgl(start_idx:end_idx, function(j) {
        groups_at_pos <- adaptive_data$Groups[[j]][[1]]
        if (is.null(groups_at_pos) || any(is.na(groups_at_pos))) return(FALSE)
        # Check if all 8 founders are in different groups (1:8)
        length(unique(groups_at_pos)) == 8 && all(sort(unique(groups_at_pos)) == 1:8)
      })
      
      # Filter to only OK positions (not all 21 positions)
      valid_haps <- window_haps[valid_positions & map_lgl(window_haps, ~ !any(is.na(.x)))]
      
      if (quality_ok && length(valid_haps) > 0) {
        # Average over only the OK positions (not all 21 positions)
        avg_haps <- reduce(valid_haps, `+`) / length(valid_haps)
        # Normalize so they sum to 1
        avg_haps <- avg_haps / sum(avg_haps)
        names(avg_haps) <- founders
        return(avg_haps)
      } else {
        # If <17 positions are OK, return NAs
        return(set_names(rep(NA, length(founders)), founders))
      }
    })),
    
    # Average the error matrices using 21-position sliding window
    Err = list(map(seq_along(sample[[1]]), function(i) {
      start_idx <- max(1, i - 10)  # 10 positions before
      end_idx <- min(nrow(adaptive_data), i + 10)  # 10 positions after
      
      # Get error matrices for this sample across the 21-position window
      window_errs <- map(start_idx:end_idx, function(j) {
        adaptive_data$Err[[j]][[1]][[i]]
      })
      
      # Only use positions where the original estimate had 8 groups
      valid_positions <- map_lgl(start_idx:end_idx, function(j) {
        groups_at_pos <- adaptive_data$Groups[[j]][[1]]
        if (is.null(groups_at_pos) || any(is.na(groups_at_pos))) return(FALSE)
        # Check if all 8 founders are in different groups (1:8)
        length(unique(groups_at_pos)) == 8 && all(sort(unique(groups_at_pos)) == 1:8)
      })
      
      # Filter to only OK positions (not all 21 positions)
      valid_errs <- window_errs[valid_positions & map_lgl(window_errs, ~ !any(is.na(.x)))]
      
      if (quality_ok && length(valid_errs) > 0) {
        # Average over only the OK positions (not all 21 positions)
        avg_err <- reduce(valid_errs, `+`) / length(valid_errs)
        rownames(avg_err) <- founders
        colnames(avg_err) <- founders
        return(avg_err)
      } else {
        # If <17 positions are OK, return NA matrix
        na_matrix <- matrix(NA, length(founders), length(founders))
        rownames(na_matrix) <- founders
        colnames(na_matrix) <- founders
        return(na_matrix)
      }
    })),
    
    # Names remain the same
    Names = list(map(seq_along(sample[[1]]), function(i) {
      return(founders)
    }))
  ) %>%
  select(-quality_count)  # Clean up temporary column

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
