#!/usr/bin/env Rscript

# Investigate what's driving the differences between h_cutoff=10 and fixed 50kb
# Are the differences due to:
# 1. Haplotype treatment differences changing rapidly?
# 2. Error matrices behaving poorly from position to position?

library(tidyverse)

# Load the reshaped data for both methods
h10_file <- "process/ZINC2_h10/adapt_h10/R.haps.chr3R.out.rds"
fixed_file <- "/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/R.haps.chr3R.out.rds"

h10_data <- readRDS(h10_file)
fixed_data <- readRDS(fixed_file)

# Load design file
design_file <- "/dfs7/adl/tdlong/fly_pool/XQTL2/helpfiles/ZINC2/Zinc2.test.M.txt"
design.df <- read.table(design_file, header = TRUE)

# Filter to test region
test_positions <- seq(20000000, 20200000, by = 10000)

h10_test <- h10_data %>% filter(pos %in% test_positions)
fixed_test <- fixed_data %>% filter(pos %in% test_positions)

cat("=== INVESTIGATING HAPLOTYPE vs ERROR CHANGES ===\n")
cat("Test region: 20,000,000-20,200,000 (21 positions)\n")
cat("Focus: Treatment differences (W vs Z) and their stability\n\n")

# Function to calculate treatment differences for a position
calculate_treatment_differences <- function(data, pos) {
  pos_data <- data %>% filter(pos == !!pos)
  
  # Extract sample data and join with design
  sample_data <- pos_data %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample == bam)) %>%
    filter(!is.na(TRT))
  
  # Calculate mean haplotype frequencies by treatment
  treatment_means <- sample_data %>%
    group_by(TRT) %>%
    summarise(
      B1_mean = mean(map_dbl(Haps, ~ .x[1])),
      B2_mean = mean(map_dbl(Haps, ~ .x[2])),
      B3_mean = mean(map_dbl(Haps, ~ .x[3])),
      B4_mean = mean(map_dbl(Haps, ~ .x[4])),
      B5_mean = mean(map_dbl(Haps, ~ .x[5])),
      B6_mean = mean(map_dbl(Haps, ~ .x[6])),
      B7_mean = mean(map_dbl(Haps, ~ .x[7])),
      AB8_mean = mean(map_dbl(Haps, ~ .x[8])),
      .groups = "drop"
    )
  
  # Calculate treatment differences (W - Z)
  if (nrow(treatment_means) == 2) {
    w_row <- treatment_means %>% filter(TRT == "W")
    z_row <- treatment_means %>% filter(TRT == "Z")
    
    if (nrow(w_row) > 0 && nrow(z_row) > 0) {
      differences <- w_row %>%
        select(starts_with("B")) - z_row %>%
        select(starts_with("B"))
      
      return(differences)
    }
  }
  
  return(NULL)
}

# Function to calculate error variance differences for a position
calculate_error_differences <- function(data, pos) {
  pos_data <- data %>% filter(pos == !!pos)
  
  # Extract sample data and join with design
  sample_data <- pos_data %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample == bam)) %>%
    filter(!is.na(TRT))
  
  # Calculate mean error variances by treatment
  treatment_errors <- sample_data %>%
    group_by(TRT) %>%
    summarise(
      B1_err = mean(map_dbl(Err, ~ diag(.x)[1])),
      B2_err = mean(map_dbl(Err, ~ diag(.x)[2])),
      B3_err = mean(map_dbl(Err, ~ diag(.x)[3])),
      B4_err = mean(map_dbl(Err, ~ diag(.x)[4])),
      B5_err = mean(map_dbl(Err, ~ diag(.x)[5])),
      B6_err = mean(map_dbl(Err, ~ diag(.x)[6])),
      B7_err = mean(map_dbl(Err, ~ diag(.x)[7])),
      AB8_err = mean(map_dbl(Err, ~ diag(.x)[8])),
      .groups = "drop"
    )
  
  # Calculate error differences (W - Z)
  if (nrow(treatment_errors) == 2) {
    w_row <- treatment_errors %>% filter(TRT == "W")
    z_row <- treatment_errors %>% filter(TRT == "Z")
    
    if (nrow(w_row) > 0 && nrow(z_row) > 0) {
      differences <- w_row %>%
        select(starts_with("B")) - z_row %>%
        select(starts_with("B"))
      
      return(differences)
    }
  }
  
  return(NULL)
}

# Calculate treatment differences for all positions
cat("=== HAPLOTYPE TREATMENT DIFFERENCES (W - Z) ===\n")
cat("H_cutoff=10:\n")

h10_hap_diffs <- map_dfr(test_positions, ~ {
  diff <- calculate_treatment_differences(h10_test, .x)
  if (!is.null(diff)) {
    diff %>% mutate(pos = .x)
  } else {
    NULL
  }
})

if (nrow(h10_hap_diffs) > 0) {
  print(h10_hap_diffs, n = Inf)
  
  # Calculate changes between adjacent positions
  h10_hap_changes <- h10_hap_diffs %>%
    arrange(pos) %>%
    mutate(
      B1_change = abs(B1_mean - lag(B1_mean)),
      B2_change = abs(B2_mean - lag(B2_mean)),
      B3_change = abs(B3_mean - lag(B3_mean)),
      B4_change = abs(B4_mean - lag(B4_mean)),
      B5_change = abs(B5_mean - lag(B5_mean)),
      B6_change = abs(B6_mean - lag(B6_mean)),
      B7_change = abs(B7_mean - lag(B7_mean)),
      AB8_change = abs(AB8_mean - lag(AB8_mean))
    ) %>%
    select(pos, starts_with("B") & ends_with("_change")) %>%
    filter(!is.na(B1_change))
  
  cat("\nChanges between adjacent positions:\n")
  print(h10_hap_changes, n = Inf)
}

cat("\nFixed 50kb:\n")

fixed_hap_diffs <- map_dfr(test_positions, ~ {
  diff <- calculate_treatment_differences(fixed_test, .x)
  if (!is.null(diff)) {
    diff %>% mutate(pos = .x)
  } else {
    NULL
  }
})

if (nrow(fixed_hap_diffs) > 0) {
  print(fixed_hap_diffs, n = Inf)
  
  # Calculate changes between adjacent positions
  fixed_hap_changes <- fixed_hap_diffs %>%
    arrange(pos) %>%
    mutate(
      B1_change = abs(B1_mean - lag(B1_mean)),
      B2_change = abs(B2_mean - lag(B2_mean)),
      B3_change = abs(B3_mean - lag(B3_mean)),
      B4_change = abs(B4_mean - lag(B4_mean)),
      B5_change = abs(B5_mean - lag(B5_mean)),
      B6_change = abs(B6_mean - lag(B6_mean)),
      B7_change = abs(B7_mean - lag(B7_mean)),
      AB8_change = abs(AB8_mean - lag(AB8_mean))
    ) %>%
    select(pos, starts_with("B") & ends_with("_change")) %>%
    filter(!is.na(B1_change))
  
  cat("\nChanges between adjacent positions:\n")
  print(fixed_hap_changes, n = Inf)
}

# Calculate error differences for all positions
cat("\n=== ERROR VARIANCE TREATMENT DIFFERENCES (W - Z) ===\n")
cat("H_cutoff=10:\n")

h10_err_diffs <- map_dfr(test_positions, ~ {
  diff <- calculate_error_differences(h10_test, .x)
  if (!is.null(diff)) {
    diff %>% mutate(pos = .x)
  } else {
    NULL
  }
})

if (nrow(h10_err_diffs) > 0) {
  print(h10_err_diffs, n = Inf)
  
  # Calculate changes between adjacent positions
  h10_err_changes <- h10_err_diffs %>%
    arrange(pos) %>%
    mutate(
      B1_err_change = abs(B1_err - lag(B1_err)),
      B2_err_change = abs(B2_err - lag(B2_err)),
      B3_err_change = abs(B3_err - lag(B3_err)),
      B4_err_change = abs(B4_err - lag(B4_err)),
      B5_err_change = abs(B5_err - lag(B5_err)),
      B6_err_change = abs(B6_err - lag(B6_err)),
      B7_err_change = abs(B7_err - lag(B7_err)),
      AB8_err_change = abs(AB8_err - lag(AB8_err))
    ) %>%
    select(pos, starts_with("B") & ends_with("_err_change")) %>%
    filter(!is.na(B1_err_change))
  
  cat("\nChanges between adjacent positions:\n")
  print(h10_err_changes, n = Inf)
}

cat("\nFixed 50kb:\n")

fixed_err_diffs <- map_dfr(test_positions, ~ {
  diff <- calculate_error_differences(fixed_test, .x)
  if (!is.null(diff)) {
    diff %>% mutate(pos = .x)
  } else {
    NULL
  }
})

if (nrow(fixed_err_diffs) > 0) {
  print(fixed_err_diffs, n = Inf)
  
  # Calculate changes between adjacent positions
  fixed_err_changes <- fixed_err_diffs %>%
    arrange(pos) %>%
    mutate(
      B1_err_change = abs(B1_err - lag(B1_err)),
      B2_err_change = abs(B2_err - lag(B2_err)),
      B3_err_change = abs(B3_err - lag(B3_err)),
      B4_err_change = abs(B4_err - lag(B4_err)),
      B5_err_change = abs(B5_err - lag(B5_err)),
      B6_err_change = abs(B6_err - lag(B6_err)),
      B7_err_change = abs(B7_err - lag(B7_err)),
      AB8_err_change = abs(AB8_err - lag(AB8_err))
    ) %>%
    select(pos, starts_with("B") & ends_with("_err_change")) %>%
    filter(!is.na(B1_err_change))
  
  cat("\nChanges between adjacent positions:\n")
  print(fixed_err_changes, n = Inf)
}

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")

if (nrow(h10_hap_changes) > 0 && nrow(fixed_hap_changes) > 0) {
  # Haplotype change summary
  h10_hap_summary <- h10_hap_changes %>%
    summarise(across(starts_with("B") & ends_with("_change"), 
                     list(mean = mean, sd = sd, max = max), na.rm = TRUE))
  
  fixed_hap_summary <- fixed_hap_changes %>%
    summarise(across(starts_with("B") & ends_with("_change"), 
                     list(mean = mean, sd = sd, max = max), na.rm = TRUE))
  
  cat("H_cutoff=10 haplotype treatment difference changes (mean ± sd, max):\n")
  print(h10_hap_summary)
  
  cat("\nFixed 50kb haplotype treatment difference changes (mean ± sd, max):\n")
  print(fixed_hap_summary)
}

if (nrow(h10_err_changes) > 0 && nrow(fixed_err_changes) > 0) {
  # Error change summary
  h10_err_summary <- h10_err_changes %>%
    summarise(across(starts_with("B") & ends_with("_err_change"), 
                     list(mean = mean, sd = sd, max = max), na.rm = TRUE))
  
  fixed_err_summary <- fixed_err_changes %>%
    summarise(across(starts_with("B") & ends_with("_err_change"), 
                     list(mean = mean, sd = sd, max = max), na.rm = TRUE))
  
  cat("\nH_cutoff=10 error variance treatment difference changes (mean ± sd, max):\n")
  print(h10_err_summary)
  
  cat("\nFixed 50kb error variance treatment difference changes (mean ± sd, max):\n")
  print(fixed_err_summary)
}