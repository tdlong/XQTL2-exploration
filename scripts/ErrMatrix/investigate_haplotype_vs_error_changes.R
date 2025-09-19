#!/usr/bin/env Rscript

# Investigate what's driving the differences between h_cutoff=10 and fixed 50kb
# Are the differences due to:
# 1. Haplotype frequencies changing rapidly?
# 2. Error matrices behaving poorly?

library(tidyverse)

# Load the reshaped data for both methods
h10_file <- "process/ZINC2_h10/adapt_h10/R.haps.chr3R.out.rds"
fixed_file <- "/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/R.haps.chr3R.out.rds"

h10_data <- readRDS(h10_file)
fixed_data <- readRDS(fixed_file)

# Filter to test region
test_positions <- seq(20000000, 20200000, by = 10000)

h10_test <- h10_data %>% filter(pos %in% test_positions)
fixed_test <- fixed_data %>% filter(pos %in% test_positions)

cat("=== INVESTIGATING HAPLOTYPE vs ERROR CHANGES ===\n")
cat("Test region: 20,000,000-20,200,000 (21 positions)\n\n")

# Function to extract haplotype frequencies and error variances for a sample
extract_sample_data <- function(data, sample_name) {
  data %>%
    mutate(
      sample_idx = map_int(sample, ~ which(.x == sample_name)),
      haps = map2(Haps, sample_idx, ~ .x[[.y]]),
      err_diag = map2(Err, sample_idx, ~ diag(.x[[.y]]))
    ) %>%
    select(pos, haps, err_diag) %>%
    unnest_wider(haps, names_sep = "_") %>%
    unnest_wider(err_diag, names_sep = "_")
}

# Extract data for a representative sample (Rep01_W_M)
sample_name <- "Rep01_W_M"

h10_sample <- extract_sample_data(h10_test, sample_name)
fixed_sample <- extract_sample_data(fixed_test, sample_name)

cat("=== HAPLOTYPE FREQUENCY CHANGES ===\n")
cat("H_cutoff=10 - Changes between adjacent positions:\n")

# Calculate changes in haplotype frequencies between adjacent positions
h10_hap_changes <- h10_sample %>%
  arrange(pos) %>%
  mutate(
    B1_change = abs(haps_B1 - lag(haps_B1)),
    B2_change = abs(haps_B2 - lag(haps_B2)),
    B3_change = abs(haps_B3 - lag(haps_B3)),
    B4_change = abs(haps_B4 - lag(haps_B4)),
    B5_change = abs(haps_B5 - lag(haps_B5)),
    B6_change = abs(haps_B6 - lag(haps_B6)),
    B7_change = abs(haps_B7 - lag(haps_B7)),
    AB8_change = abs(haps_AB8 - lag(haps_AB8))
  ) %>%
  select(pos, starts_with("B") & ends_with("_change")) %>%
  filter(!is.na(B1_change))

print(h10_hap_changes, n = Inf)

cat("\nFixed 50kb - Changes between adjacent positions:\n")

fixed_hap_changes <- fixed_sample %>%
  arrange(pos) %>%
  mutate(
    B1_change = abs(haps_B1 - lag(haps_B1)),
    B2_change = abs(haps_B2 - lag(haps_B2)),
    B3_change = abs(haps_B3 - lag(haps_B3)),
    B4_change = abs(haps_B4 - lag(haps_B4)),
    B5_change = abs(haps_B5 - lag(haps_B5)),
    B6_change = abs(haps_B6 - lag(haps_B6)),
    B7_change = abs(haps_B7 - lag(haps_B7)),
    AB8_change = abs(haps_AB8 - lag(haps_AB8))
  ) %>%
  select(pos, starts_with("B") & ends_with("_change")) %>%
  filter(!is.na(B1_change))

print(fixed_hap_changes, n = Inf)

cat("\n=== ERROR VARIANCE CHANGES ===\n")
cat("H_cutoff=10 - Changes in error variances between adjacent positions:\n")

h10_err_changes <- h10_sample %>%
  arrange(pos) %>%
  mutate(
    B1_err_change = abs(err_diag_B1 - lag(err_diag_B1)),
    B2_err_change = abs(err_diag_B2 - lag(err_diag_B2)),
    B3_err_change = abs(err_diag_B3 - lag(err_diag_B3)),
    B4_err_change = abs(err_diag_B4 - lag(err_diag_B4)),
    B5_err_change = abs(err_diag_B5 - lag(err_diag_B5)),
    B6_err_change = abs(err_diag_B6 - lag(err_diag_B6)),
    B7_err_change = abs(err_diag_B7 - lag(err_diag_B7)),
    AB8_err_change = abs(err_diag_AB8 - lag(err_diag_AB8))
  ) %>%
  select(pos, starts_with("B") & ends_with("_err_change")) %>%
  filter(!is.na(B1_err_change))

print(h10_err_changes, n = Inf)

cat("\nFixed 50kb - Changes in error variances between adjacent positions:\n")

fixed_err_changes <- fixed_sample %>%
  arrange(pos) %>%
  mutate(
    B1_err_change = abs(err_diag_B1 - lag(err_diag_B1)),
    B2_err_change = abs(err_diag_B2 - lag(err_diag_B2)),
    B3_err_change = abs(err_diag_B3 - lag(err_diag_B3)),
    B4_err_change = abs(err_diag_B4 - lag(err_diag_B4)),
    B5_err_change = abs(err_diag_B5 - lag(err_diag_B5)),
    B6_err_change = abs(err_diag_B6 - lag(err_diag_B6)),
    B7_err_change = abs(err_diag_B7 - lag(err_diag_B7)),
    AB8_err_change = abs(err_diag_AB8 - lag(err_diag_AB8))
  ) %>%
  select(pos, starts_with("B") & ends_with("_err_change")) %>%
  filter(!is.na(B1_err_change))

print(fixed_err_changes, n = Inf)

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")

# Haplotype change summary
h10_hap_summary <- h10_hap_changes %>%
  summarise(across(starts_with("B") & ends_with("_change"), 
                   list(mean = mean, sd = sd, max = max), na.rm = TRUE))

fixed_hap_summary <- fixed_hap_changes %>%
  summarise(across(starts_with("B") & ends_with("_change"), 
                   list(mean = mean, sd = sd, max = max), na.rm = TRUE))

cat("H_cutoff=10 haplotype changes (mean ± sd, max):\n")
print(h10_hap_summary)

cat("\nFixed 50kb haplotype changes (mean ± sd, max):\n")
print(fixed_hap_summary)

# Error change summary
h10_err_summary <- h10_err_changes %>%
  summarise(across(starts_with("B") & ends_with("_err_change"), 
                   list(mean = mean, sd = sd, max = max), na.rm = TRUE))

fixed_err_summary <- fixed_err_changes %>%
  summarise(across(starts_with("B") & ends_with("_err_change"), 
                   list(mean = mean, sd = sd, max = max), na.rm = TRUE))

cat("\nH_cutoff=10 error variance changes (mean ± sd, max):\n")
print(h10_err_summary)

cat("\nFixed 50kb error variance changes (mean ± sd, max):\n")
print(fixed_err_summary)
