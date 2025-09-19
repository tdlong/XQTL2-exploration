#!/usr/bin/env Rscript

# Script to investigate smoothness of haplotype estimates across positions
# Focus on chr3R positions 20,000,000 to 20,200,000 (20 positions, 10kb apart)

library(tidyverse)

# Load the adaptive results for chr3R
adaptive_file <- "process/ZINC2_h10/adaptive_window_h10_results_chr3R.RDS"
adaptive_data <- readRDS(adaptive_file)

# Define the test region
test_positions <- seq(20000000, 20200000, by = 10000)  # 20 positions
test_sample <- "Rep01_W_F"  # Pick an arbitrary sample

# Filter to our test region and sample
test_data <- adaptive_data %>%
  filter(pos %in% test_positions, sample == test_sample) %>%
  arrange(pos)

cat("Found", nrow(test_data), "positions for sample", test_sample, "\n")
cat("Positions:", paste(test_data$pos, collapse = ", "), "\n\n")

# Extract haplotype frequencies
cat("=== HAPLOTYPE FREQUENCIES ===\n")
hap_freqs <- test_data %>%
  select(pos, Haps) %>%
  mutate(
    B1 = map_dbl(Haps, ~ .x["B1"]),
    B2 = map_dbl(Haps, ~ .x["B2"]),
    B3 = map_dbl(Haps, ~ .x["B3"]),
    B4 = map_dbl(Haps, ~ .x["B4"]),
    B5 = map_dbl(Haps, ~ .x["B5"]),
    B6 = map_dbl(Haps, ~ .x["B6"]),
    B7 = map_dbl(Haps, ~ .x["B7"]),
    B8 = map_dbl(Haps, ~ .x["B8"])
  ) %>%
  select(pos, B1, B2, B3, B4, B5, B6, B7, B8)

print(hap_freqs)

# Extract error variances (diagonal of error matrix)
cat("\n=== ERROR VARIANCES (sqrt of diagonal) ===\n")
error_vars <- test_data %>%
  select(pos, Err) %>%
  mutate(
    sqrt_var_B1 = map_dbl(Err, ~ sqrt(.x["B1", "B1"])),
    sqrt_var_B2 = map_dbl(Err, ~ sqrt(.x["B2", "B2"])),
    sqrt_var_B3 = map_dbl(Err, ~ sqrt(.x["B3", "B3"])),
    sqrt_var_B4 = map_dbl(Err, ~ sqrt(.x["B4", "B4"])),
    sqrt_var_B5 = map_dbl(Err, ~ sqrt(.x["B5", "B5"])),
    sqrt_var_B6 = map_dbl(Err, ~ sqrt(.x["B6", "B6"])),
    sqrt_var_B7 = map_dbl(Err, ~ sqrt(.x["B7", "B7"])),
    sqrt_var_B8 = map_dbl(Err, ~ sqrt(.x["B8", "B8"]))
  ) %>%
  select(pos, sqrt_var_B1, sqrt_var_B2, sqrt_var_B3, sqrt_var_B4, 
         sqrt_var_B5, sqrt_var_B6, sqrt_var_B7, sqrt_var_B8)

print(error_vars)

# Calculate differences between adjacent positions
cat("\n=== HAPLOTYPE FREQUENCY DIFFERENCES (adjacent positions) ===\n")
hap_diffs <- hap_freqs %>%
  mutate(
    across(B1:B8, ~ .x - lag(.x), .names = "diff_{.col}")
  ) %>%
  select(pos, starts_with("diff_")) %>%
  filter(!is.na(diff_B1))

print(hap_diffs)

# Calculate differences in error variances
cat("\n=== ERROR VARIANCE DIFFERENCES (adjacent positions) ===\n")
err_diffs <- error_vars %>%
  mutate(
    across(sqrt_var_B1:sqrt_var_B8, ~ .x - lag(.x), .names = "diff_{.col}")
  ) %>%
  select(pos, starts_with("diff_")) %>%
  filter(!is.na(diff_sqrt_var_B1))

print(err_diffs)

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Max haplotype frequency change:", max(abs(hap_diffs[,2:9]), na.rm = TRUE), "\n")
cat("Mean haplotype frequency change:", mean(abs(hap_diffs[,2:9]), na.rm = TRUE), "\n")
cat("Max error variance change:", max(abs(err_diffs[,2:9]), na.rm = TRUE), "\n")
cat("Mean error variance change:", mean(abs(err_diffs[,2:9]), na.rm = TRUE), "\n")

# Additional analysis: check for patterns
cat("\n=== PATTERN ANALYSIS ===\n")
cat("Number of positions with large haplotype changes (>0.1):", 
    sum(abs(hap_diffs[,2:9]) > 0.1, na.rm = TRUE), "\n")
cat("Number of positions with large error changes (>0.01):", 
    sum(abs(err_diffs[,2:9]) > 0.01, na.rm = TRUE), "\n")

# Check if there are any NA values
cat("Positions with NA haplotype frequencies:", sum(is.na(hap_freqs[,2:9])), "\n")
cat("Positions with NA error variances:", sum(is.na(error_vars[,2:9])), "\n")
