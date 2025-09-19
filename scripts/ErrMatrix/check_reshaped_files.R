#!/usr/bin/env Rscript

# Script to check reshaped files for consistency
# 1. Do reshaped files give same results as long format?
# 2. Are samples in same order at every position?

library(tidyverse)

# Load both long and reshaped files
long_file <- "process/ZINC2_h10/adaptive_window_h10_results_chr3R.RDS"
reshaped_file <- "process/ZINC2_h10/adapt_h10/R.haps.chr3R.out.rds"

long_data <- readRDS(long_file)
reshaped_data <- readRDS(reshaped_file)

cat("=== CHECKING RESHAPED FILE CONSISTENCY ===\n")
cat("Long format dimensions:", dim(long_data), "\n")
cat("Reshaped format dimensions:", dim(reshaped_data), "\n\n")

# Define test region
test_positions <- seq(20000000, 20200000, by = 10000)
test_sample <- "Rep01_W_M"  # Use male sample

# 1. Check if reshaped files give same results
cat("=== 1. COMPARING LONG vs RESHAPED FORMAT ===\n")

# Extract from long format
long_subset <- long_data %>%
  filter(pos %in% test_positions, sample == test_sample) %>%
  arrange(pos) %>%
  select(pos, Haps) %>%
  mutate(
    B1 = map_dbl(Haps, ~ .x["B1"]),
    B2 = map_dbl(Haps, ~ .x["B2"]),
    B3 = map_dbl(Haps, ~ .x["B3"]),
    B4 = map_dbl(Haps, ~ .x["B4"]),
    B5 = map_dbl(Haps, ~ .x["B5"]),
    B6 = map_dbl(Haps, ~ .x["B6"]),
    B7 = map_dbl(Haps, ~ .x["B7"]),
    AB8 = map_dbl(Haps, ~ .x["AB8"])
  ) %>%
  select(pos, B1, B2, B3, B4, B5, B6, B7, AB8)

# Extract from reshaped format
reshaped_subset <- reshaped_data %>%
  filter(pos %in% test_positions) %>%
  arrange(pos) %>%
  select(pos, sample, Haps) %>%
  # Find the sample index for our test sample
  mutate(sample_index = map_int(sample, ~ which(.x == test_sample))) %>%
  filter(!is.na(sample_index)) %>%
  mutate(
    B1 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B1"]),
    B2 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B2"]),
    B3 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B3"]),
    B4 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B4"]),
    B5 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B5"]),
    B6 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B6"]),
    B7 = map_dbl(Haps, ~ .x[[sample_index[1]]]["B7"]),
    AB8 = map_dbl(Haps, ~ .x[[sample_index[1]]]["AB8"])
  ) %>%
  select(pos, B1, B2, B3, B4, B5, B6, B7, AB8)

# Compare the two
comparison <- long_subset %>%
  left_join(reshaped_subset, by = "pos", suffix = c("_long", "_reshaped")) %>%
  mutate(
    diff_B1 = abs(B1_long - B1_reshaped),
    diff_B2 = abs(B2_long - B2_reshaped),
    diff_B3 = abs(B3_long - B3_reshaped),
    diff_B4 = abs(B4_long - B4_reshaped),
    diff_B5 = abs(B5_long - B5_reshaped),
    diff_B6 = abs(B6_long - B6_reshaped),
    diff_B7 = abs(B7_long - B7_reshaped),
    diff_AB8 = abs(AB8_long - AB8_reshaped)
  )

cat("Maximum differences between long and reshaped formats:\n")
max_diffs <- comparison %>%
  summarise(
    B1 = max(diff_B1, na.rm = TRUE),
    B2 = max(diff_B2, na.rm = TRUE),
    B3 = max(diff_B3, na.rm = TRUE),
    B4 = max(diff_B4, na.rm = TRUE),
    B5 = max(diff_B5, na.rm = TRUE),
    B6 = max(diff_B6, na.rm = TRUE),
    B7 = max(diff_B7, na.rm = TRUE),
    B8 = max(diff_AB8, na.rm = TRUE)
  )
print(max_diffs)

# 2. Check sample order consistency
cat("\n=== 2. CHECKING SAMPLE ORDER CONSISTENCY ===\n")

# Get sample order from first few positions
sample_orders <- reshaped_data %>%
  filter(pos %in% test_positions[1:5]) %>%  # Check first 5 positions
  select(pos, sample) %>%
  mutate(sample_list = map(sample, ~ .x)) %>%
  select(pos, sample_list)

cat("Sample order at first 5 positions:\n")
for(i in 1:nrow(sample_orders)) {
  cat("Position", sample_orders$pos[i], ":", paste(sample_orders$sample_list[[i]], collapse = " "), "\n")
}

# Check if sample order is identical across positions
first_order <- sample_orders$sample_list[[1]]
all_identical <- all(map_lgl(sample_orders$sample_list, ~ identical(.x, first_order)))

cat("\nAre all sample orders identical?", all_identical, "\n")

if(!all_identical) {
  cat("WARNING: Sample order varies across positions!\n")
  cat("This could cause issues if your software assumes constant order.\n")
} else {
  cat("Sample order is consistent across all positions.\n")
}
