#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

# =============================================================================
# Create Smooth Haplotype Estimator
# =============================================================================
# 
# This script creates a "smooth_h4" estimator by applying 21-position sliding
# window smoothing to adaptive_h4 results, then rescaling to ensure founder
# frequencies sum to 1.
#
# USAGE:
# Rscript scripts/production/create_smooth_haplotype_estimator.R <chr> <param_file> <output_dir>
#
# EXAMPLE:
# Rscript scripts/production/create_smooth_haplotype_estimator.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/production/create_smooth_haplotype_estimator.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

cat("=== CREATING SMOOTH HAPLOTYPE ESTIMATOR ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

# Load the adaptive_h4 results
adaptive_file <- file.path(results_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))

if (!file.exists(adaptive_file)) {
  stop("Adaptive h4 results not found: ", adaptive_file)
}

cat("Loading adaptive_h4 results...\n")
adaptive_data <- readRDS(adaptive_file)

cat("✓ Adaptive h4 data loaded:", nrow(adaptive_data), "rows\n")
cat("Samples available:", paste(unique(adaptive_data$sample), collapse = ", "), "\n")
cat("Founders available:", paste(unique(adaptive_data$founder), collapse = ", "), "\n\n")

# Function to apply sliding window smoothing with quality control
smooth_founder_frequencies <- function(data, founder_name) {
  cat("Smoothing founder:", founder_name, "\n")
  
  # Filter to this founder and arrange by position
  founder_data <- data %>%
    filter(founder == founder_name) %>%
    arrange(pos)
  
  if (nrow(founder_data) == 0) {
    cat("  No data for founder:", founder_name, "\n")
    return(founder_data)
  }
  
  # Apply sliding window smoothing
  founder_data <- founder_data %>%
    group_by(sample) %>%
    arrange(pos) %>%
    mutate(
      # Smooth the frequency values (only over estimate_OK==1 positions)
      smoothed_freq = zoo::rollapply(
        data = frequency,
        width = 21,
        FUN = function(x) {
          # Only use positions where estimate_OK==1
          ok_positions <- which(!is.na(x) & founder_data$estimate_OK[seq_along(x)] == 1)
          if (length(ok_positions) == 0) return(NA)
          mean(x[ok_positions], na.rm = TRUE)
        },
        fill = NA,
        align = "center",
        partial = FALSE
      ),
      # Calculate quality: how many positions in window are OK?
      quality_count = zoo::rollapply(
        data = estimate_OK,
        width = 21,
        FUN = function(x) sum(x == 1, na.rm = TRUE),
        fill = NA,
        align = "center",
        partial = FALSE
      ),
      # New estimate_OK: OK if at least 17 out of 21 positions are OK
      new_estimate_OK = ifelse(quality_count >= 17, 1, 0)
    ) %>%
    ungroup()
  
  cat("  Smoothed", nrow(founder_data), "positions\n")
  cat("  Quality: ", sum(founder_data$new_estimate_OK == 1, na.rm = TRUE), "OK, ", 
      sum(founder_data$new_estimate_OK == 0, na.rm = TRUE), "not OK\n")
  
  return(founder_data)
}

# Apply smoothing to each founder
cat("Applying 21-position sliding window smoothing to each founder...\n")
smoothed_data <- adaptive_data %>%
  group_by(founder) %>%
  group_modify(~ smooth_founder_frequencies(.x, .y$founder[1])) %>%
  ungroup()

# Now rescale each position so all founder frequencies sum to 1
cat("\nRescaling founder frequencies to sum to 1...\n")
rescaled_data <- smoothed_data %>%
  group_by(sample, pos) %>%
  mutate(
    # Calculate sum of all founder frequencies at this position
    total_freq = sum(smoothed_freq, na.rm = TRUE),
    # Rescale each founder frequency
    rescaled_freq = ifelse(total_freq > 0, smoothed_freq / total_freq, NA),
    # Keep the new estimate_OK
    estimate_OK = new_estimate_OK
  ) %>%
  ungroup() %>%
  # Clean up intermediate columns
  select(-smoothed_freq, -quality_count, -new_estimate_OK, -total_freq) %>%
  # Rename the rescaled frequency back to frequency
  rename(frequency = rescaled_freq)

# Add method information
rescaled_data <- rescaled_data %>%
  mutate(method = "smooth_h4")

# Save the smooth_h4 results
output_file <- file.path(results_dir, paste0("smooth_h4_results_", chr, ".RDS"))
saveRDS(rescaled_data, output_file)

cat("✓ Smooth h4 results saved to:", output_file, "\n")

# Print summary statistics
cat("\n=== SMOOTH H4 SUMMARY ===\n")
cat("Total positions:", nrow(rescaled_data), "\n")
cat("Positions with estimate_OK = 1:", sum(rescaled_data$estimate_OK == 1, na.rm = TRUE), "\n")
cat("Positions with estimate_OK = 0:", sum(rescaled_data$estimate_OK == 0, na.rm = TRUE), "\n")
cat("Positions with estimate_OK = NA:", sum(is.na(rescaled_data$estimate_OK)), "\n")

# Check that frequencies sum to 1
freq_sums <- rescaled_data %>%
  group_by(sample, pos) %>%
  summarize(total_freq = sum(frequency, na.rm = TRUE), .groups = "drop")

cat("\nFrequency sum validation:\n")
cat("Mean sum:", round(mean(freq_sums$total_freq, na.rm = TRUE), 6), "\n")
cat("Min sum:", round(min(freq_sums$total_freq, na.rm = TRUE), 6), "\n")
cat("Max sum:", round(max(freq_sums$total_freq, na.rm = TRUE), 6), "\n")

cat("\n✓ Smooth h4 estimator creation complete!\n")
