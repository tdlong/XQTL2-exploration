#!/usr/bin/env Rscript

# Script to investigate smoothness of haplotype estimates across positions
# Focus on chr3R positions 20,000,000 to 20,200,000 (20 positions, 10kb apart)

library(tidyverse)

# Load both h10 and h4 results for comparison
adaptive_file_h10 <- "process/ZINC2_h10/adaptive_window_h10_results_chr3R.RDS"
adaptive_file_h4 <- "process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_chr3R.RDS"

adaptive_data_h10 <- readRDS(adaptive_file_h10)
adaptive_data_h4 <- readRDS(adaptive_file_h4)

# Define the test region
test_positions <- seq(20000000, 20200000, by = 10000)  # 20 positions
test_sample <- "Rep01_W_F"  # Pick an arbitrary sample

# Filter to our test region and sample for both datasets
test_data_h10 <- adaptive_data_h10 %>%
  filter(pos %in% test_positions, sample == test_sample) %>%
  arrange(pos)

test_data_h4 <- adaptive_data_h4 %>%
  filter(pos %in% test_positions, sample == test_sample) %>%
  arrange(pos)

cat("Found", nrow(test_data_h10), "positions for sample", test_sample, "in h10 data\n")
cat("Found", nrow(test_data_h4), "positions for sample", test_sample, "in h4 data\n")
cat("Positions:", paste(test_data_h10$pos, collapse = ", "), "\n\n")

# Function to extract and format data
extract_data <- function(test_data, label) {
  first_hap <- test_data$Haps[[1]]
  founder_names <- names(first_hap)
  
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
      AB8 = map_dbl(Haps, ~ .x["AB8"])
    ) %>%
    select(pos, B1, B2, B3, B4, B5, B6, B7, AB8) %>%
    mutate(across(B1:AB8, ~ round(.x * 1000, 0)))
  
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
      AB8 = map_dbl(Err, ~ sqrt(.x["AB8", "AB8"]))
    ) %>%
    select(pos, B1, B2, B3, B4, B5, B6, B7, AB8) %>%
    mutate(across(B1:AB8, ~ round(.x * 1000, 0)))
  
  return(list(hap_freqs = hap_freqs, error_vars = error_vars))
}

# Extract data for both h10 and h4
data_h10 <- extract_data(test_data_h10, "h10")
data_h4 <- extract_data(test_data_h4, "h4")

# Show h10 data
cat("=== H10 HAPLOTYPE FREQUENCIES (×1000) ===\n")
print(data_h10$hap_freqs, n = Inf)

cat("\n=== H10 ERROR VARIANCES (sqrt × 1000) ===\n")
print(data_h10$error_vars, n = Inf)

# Show h4 data
cat("\n=== H4 HAPLOTYPE FREQUENCIES (×1000) ===\n")
print(data_h4$hap_freqs, n = Inf)

cat("\n=== H4 ERROR VARIANCES (sqrt × 1000) ===\n")
print(data_h4$error_vars, n = Inf)

