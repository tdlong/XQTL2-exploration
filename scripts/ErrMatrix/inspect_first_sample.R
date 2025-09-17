#!/usr/bin/env Rscript

# Inspect raw haplotype estimates for first position and first sample
# to understand the magnitude of differences

library(tidyverse)

# Load the data
data_file <- "testing_positions_comparison.rds"
xx <- readRDS(data_file)

cat("Data structure:\n")
str(xx)

# Get first position
first_pos <- min(xx$pos)
cat("\nFirst position:", first_pos, "\n")

# Filter to first position
pos_data <- xx %>%
  dplyr::filter(pos == first_pos)

cat("\nMethods for first position:", paste(unique(pos_data$method), collapse = ", "), "\n")

# Extract first sample from each method
adapt_data <- pos_data %>%
  dplyr::filter(method == "adapt") %>%
  dplyr::pull(sample) %>%
  .[[1]]

fixed_data <- pos_data %>%
  dplyr::filter(method == "fixed") %>%
  dplyr::pull(sample) %>%
  .[[1]]

cat("\nFirst sample (adapt):", adapt_data[1], "\n")
cat("First sample (fixed):", fixed_data[1], "\n")

# Extract haplotype estimates for first sample
adapt_hap <- pos_data %>%
  dplyr::filter(method == "adapt") %>%
  dplyr::pull(Haps) %>%
  .[[1]] %>%
  .[[1]]  # First sample's haplotypes

fixed_hap <- pos_data %>%
  dplyr::filter(method == "fixed") %>%
  dplyr::pull(Haps) %>%
  .[[1]] %>%
  .[[1]]  # First sample's haplotypes

cat("\nHaplotype estimates for first sample at position", first_pos, ":\n")
cat("Adapt:", paste(round(adapt_hap, 6), collapse = ", "), "\n")
cat("Fixed:", paste(round(fixed_hap, 6), collapse = ", "), "\n")

# Calculate differences
hap_diff <- adapt_hap - fixed_hap
cat("\nDifferences (adapt - fixed):", paste(round(hap_diff, 6), collapse = ", "), "\n")
cat("Sum of absolute differences:", round(sum(abs(hap_diff)), 6), "\n")

# Show as data frame for easier reading
comparison_df <- tibble(
  haplotype = 1:8,
  adapt = adapt_hap,
  fixed = fixed_hap,
  difference = hap_diff,
  abs_difference = abs(hap_diff)
)

cat("\nDetailed comparison:\n")
print(comparison_df)

# Also show error matrix diagonals for comparison
adapt_err_diag <- pos_data %>%
  dplyr::filter(method == "adapt") %>%
  dplyr::pull(Err) %>%
  .[[1]] %>%
  .[[1]] %>%  # First sample's error matrix
  diag()

fixed_err_diag <- pos_data %>%
  dplyr::filter(method == "fixed") %>%
  dplyr::pull(Err) %>%
  .[[1]] %>%
  .[[1]] %>%  # First sample's error matrix
  diag()

cat("\nError matrix diagonals for first sample at position", first_pos, ":\n")
cat("Adapt:", paste(round(adapt_err_diag, 6), collapse = ", "), "\n")
cat("Fixed:", paste(round(fixed_err_diag, 6), collapse = ", "), "\n")
cat("Sum adapt:", round(sum(adapt_err_diag), 6), "\n")
cat("Sum fixed:", round(sum(fixed_err_diag), 6), "\n")
cat("Difference (adapt - fixed):", round(sum(adapt_err_diag) - sum(fixed_err_diag), 6), "\n")
