#!/usr/bin/env Rscript

# Script to analyze treatment differences in haplotype frequencies
# Calculate average difference (W - Z) for each founder at each position
# Separate analysis for males and females

library(tidyverse)

# Function to extract and calculate treatment differences
analyze_treatment_differences <- function(adaptive_file, dataset_name) {
  cat("=== ANALYZING", dataset_name, "===\\n")
  
  # Load data
  adaptive_data <- readRDS(adaptive_file)
  
  # Define the test region
  test_positions <- seq(20000000, 20200000, by = 10000)
  
  # Filter to our test region
  test_data <- adaptive_data %>%
    filter(pos %in% test_positions) %>%
    arrange(pos, sample)
  
  cat("Found", nrow(test_data), "observations in test region\\n")
  cat("Samples:", paste(unique(test_data$sample), collapse = ", "), "\\n\\n")
  
  # Extract haplotype frequencies
  hap_data <- test_data %>%
    select(pos, sample, Haps) %>%
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
    select(pos, sample, B1, B2, B3, B4, B5, B6, B7, AB8)
  
  # Parse sample names to get treatment and sex
  hap_data <- hap_data %>%
    mutate(
      treatment = ifelse(str_detect(sample, "_W_"), "W", "Z"),
      sex = ifelse(str_detect(sample, "_M"), "M", "F"),
      replicate = str_extract(sample, "Rep\\d+")
    )
  
  # Calculate treatment differences for males
  male_diffs <- hap_data %>%
    filter(sex == "M") %>%
    select(pos, replicate, treatment, B1:AB8) %>%
    pivot_wider(names_from = treatment, values_from = B1:AB8) %>%
    mutate(
      diff_B1 = B1_W - B1_Z,
      diff_B2 = B2_W - B2_Z,
      diff_B3 = B3_W - B3_Z,
      diff_B4 = B4_W - B4_Z,
      diff_B5 = B5_W - B5_Z,
      diff_B6 = B6_W - B6_Z,
      diff_B7 = B7_W - B7_Z,
      diff_AB8 = AB8_W - AB8_Z
    ) %>%
    select(pos, replicate, starts_with("diff_"))
  
  # Calculate average differences across replicates for each position
  male_avg_diffs <- male_diffs %>%
    group_by(pos) %>%
    summarise(
      avg_diff_B1 = mean(diff_B1, na.rm = TRUE),
      avg_diff_B2 = mean(diff_B2, na.rm = TRUE),
      avg_diff_B3 = mean(diff_B3, na.rm = TRUE),
      avg_diff_B4 = mean(diff_B4, na.rm = TRUE),
      avg_diff_B5 = mean(diff_B5, na.rm = TRUE),
      avg_diff_B6 = mean(diff_B6, na.rm = TRUE),
      avg_diff_B7 = mean(diff_B7, na.rm = TRUE),
      avg_diff_AB8 = mean(diff_AB8, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(starts_with("avg_diff_"), ~ round(.x * 1000, 0)))
  
  # Calculate treatment differences for females
  female_diffs <- hap_data %>%
    filter(sex == "F") %>%
    select(pos, replicate, treatment, B1:AB8) %>%
    pivot_wider(names_from = treatment, values_from = B1:AB8) %>%
    mutate(
      diff_B1 = B1_W - B1_Z,
      diff_B2 = B2_W - B2_Z,
      diff_B3 = B3_W - B3_Z,
      diff_B4 = B4_W - B4_Z,
      diff_B5 = B5_W - B5_Z,
      diff_B6 = B6_W - B6_Z,
      diff_B7 = B7_W - B7_Z,
      diff_AB8 = AB8_W - AB8_Z
    ) %>%
    select(pos, replicate, starts_with("diff_"))
  
  # Calculate average differences across replicates for each position
  female_avg_diffs <- female_diffs %>%
    group_by(pos) %>%
    summarise(
      avg_diff_B1 = mean(diff_B1, na.rm = TRUE),
      avg_diff_B2 = mean(diff_B2, na.rm = TRUE),
      avg_diff_B3 = mean(diff_B3, na.rm = TRUE),
      avg_diff_B4 = mean(diff_B4, na.rm = TRUE),
      avg_diff_B5 = mean(diff_B5, na.rm = TRUE),
      avg_diff_B6 = mean(diff_B6, na.rm = TRUE),
      avg_diff_B7 = mean(diff_B7, na.rm = TRUE),
      avg_diff_AB8 = mean(diff_AB8, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(starts_with("avg_diff_"), ~ round(.x * 1000, 0)))
  
  return(list(male = male_avg_diffs, female = female_avg_diffs))
}

# Analyze both datasets
h10_results <- analyze_treatment_differences(
  "process/ZINC2_h10/adaptive_window_h10_results_chr3R.RDS", 
  "H10"
)

h4_results <- analyze_treatment_differences(
  "process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_chr3R.RDS", 
  "H4"
)

# Display results
cat("\\n=== H10 MALE TREATMENT DIFFERENCES (W - Z, ×1000) ===\\n")
print(h10_results$male, n = Inf)

cat("\\n=== H10 FEMALE TREATMENT DIFFERENCES (W - Z, ×1000) ===\\n")
print(h10_results$female, n = Inf)

cat("\\n=== H4 MALE TREATMENT DIFFERENCES (W - Z, ×1000) ===\\n")
print(h4_results$male, n = Inf)

cat("\\n=== H4 FEMALE TREATMENT DIFFERENCES (W - Z, ×1000) ===\\n")
print(h4_results$female, n = Inf)
