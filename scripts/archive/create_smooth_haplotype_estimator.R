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
cat("Founder columns:", paste(c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8"), collapse = ", "), "\n\n")

# Get founder column names
founder_cols <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Apply sliding window smoothing to each founder column
cat("Applying 21-position sliding window smoothing to each founder...\n")
smoothed_data <- adaptive_data %>%
  group_by(sample) %>%
  arrange(pos) %>%
  mutate(
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

# Apply smoothing to each founder column
for (founder in founder_cols) {
  cat("Smoothing founder:", founder, "\n")
  
  smoothed_data <- smoothed_data %>%
    group_by(sample) %>%
    arrange(pos) %>%
    mutate(
      # Create smoothed version of this founder
      !!paste0(founder, "_smoothed") := zoo::rollapply(
        data = !!sym(founder),
        width = 21,
        FUN = function(x) {
          # Only use positions where estimate_OK==1
          ok_positions <- which(!is.na(x) & smoothed_data$estimate_OK[seq_along(x)] == 1)
          if (length(ok_positions) == 0) return(NA)
          mean(x[ok_positions], na.rm = TRUE)
        },
        fill = NA,
        align = "center",
        partial = FALSE
      )
    ) %>%
    ungroup()
}

# Now rescale each position so all founder frequencies sum to 1
cat("\nRescaling founder frequencies to sum to 1...\n")

# Calculate total frequency for each position (sum of all smoothed founder frequencies)
rescaled_data <- smoothed_data %>%
  mutate(
    total_freq = B1_smoothed + B2_smoothed + B3_smoothed + B4_smoothed + 
                 B5_smoothed + B6_smoothed + B7_smoothed + AB8_smoothed
  )

# Rescale each founder frequency
for (founder in founder_cols) {
  smoothed_col <- paste0(founder, "_smoothed")
  rescaled_data[[founder]] <- ifelse(
    rescaled_data$total_freq > 0, 
    rescaled_data[[smoothed_col]] / rescaled_data$total_freq, 
    NA
  )
}

# Clean up and finalize
rescaled_data <- rescaled_data %>%
  select(-ends_with("_smoothed"), -quality_count, -total_freq) %>%
  mutate(
    estimate_OK = new_estimate_OK,
    method = "smooth_h4"
  ) %>%
  select(-new_estimate_OK)

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
  mutate(total_freq = B1 + B2 + B3 + B4 + B5 + B6 + B7 + AB8) %>%
  select(sample, pos, total_freq)

cat("\nFrequency sum validation:\n")
cat("Mean sum:", round(mean(freq_sums$total_freq, na.rm = TRUE), 6), "\n")
cat("Min sum:", round(min(freq_sums$total_freq, na.rm = TRUE), 6), "\n")
cat("Max sum:", round(max(freq_sums$total_freq, na.rm = TRUE), 6), "\n")

cat("\n✓ Smooth h4 estimator creation complete!\n")
