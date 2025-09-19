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

# Extract haplotype frequencies - use actual founder names
cat("=== HAPLOTYPE FREQUENCIES (×1000) ===\n")
first_hap <- test_data$Haps[[1]]
founder_names <- names(first_hap)
cat("Founders:", paste(founder_names, collapse = " "), "\n")

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
    AB8 = map_dbl(Haps, ~ .x["AB8"])  # Use AB8, not B8
  ) %>%
  select(pos, B1, B2, B3, B4, B5, B6, B7, AB8) %>%
  mutate(across(B1:AB8, ~ round(.x * 1000, 0)))

print(hap_freqs, n = Inf)

# Extract error variances (sqrt of diagonal) - use actual founder names
cat("\n=== ERROR VARIANCES (sqrt × 1000) ===\n")
error_vars <- test_data %>%
  select(pos, Err) %>%
  mutate(
    B1 = map_dbl(Err, ~ sqrt(.x["B1", "B1"])),
    B2 = map_dbl(Err, ~ sqrt(.x["B2", "B2"])),
    B3 = map_dbl(Err, ~ sqrt(.x["B3", "B3"])),
    B4 = map_dbl(Err, ~ sqrt(.x["B4", "B4"])),
    B5 = map_dbl(Err, ~ sqrt(.x["B5", "B5"])),
    B6 = map_dbl(Err, ~ sqrt(.x["B6", "B6"])),
    B7 = map_dbl(Err, ~ sqrt(.x["B7", "B7"])),
    AB8 = map_dbl(Err, ~ sqrt(.x["AB8", "AB8"]))  # Use AB8, not B8
  ) %>%
  select(pos, B1, B2, B3, B4, B5, B6, B7, AB8) %>%
  mutate(across(B1:AB8, ~ round(.x * 1000, 0)))

print(error_vars, n = Inf)

