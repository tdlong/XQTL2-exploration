#!/usr/bin/env Rscript

# Compare hunk data results vs production results
# Run locally after downloading both files

suppressPackageStartupMessages({
  library(tidyverse)
})

# Load the estimator function
source("scripts/ErrMatrix/archive/BASE_VAR_WIDE.R")

# Load hunk data (from extract_hunk.r)
hunk_file <- "hunk_data_chr3R_19780000_Rep01_W_F_h4.rds"
if (!file.exists(hunk_file)) {
  stop("Hunk data file not found: ", hunk_file)
}

hunk_data <- readRDS(hunk_file)
cat("=== HUNK DATA ===\n")
cat("df3 dimensions:", nrow(hunk_data$df3), "x", ncol(hunk_data$df3), "\n")
cat("Args:", paste(names(hunk_data$args), collapse = ", "), "\n\n")

# Run estimator on hunk data
cat("=== RUNNING ESTIMATOR ON HUNK DATA ===\n")
hunk_result <- do.call(est_haps_var, c(list(df3 = hunk_data$df3), hunk_data$args))

cat("Hunk result structure:\n")
str(hunk_result)

# Load production data (from testing_positions_comparison.rds)
prod_file <- "testing_positions_comparison.rds"
if (!file.exists(prod_file)) {
  stop("Production data file not found: ", prod_file)
}

prod_data <- readRDS(prod_file)
cat("\n=== PRODUCTION DATA ===\n")
cat("Total rows:", nrow(prod_data), "\n")

# Find the specific position and sample in production data
prod_row <- prod_data %>%
  dplyr::filter(pos == 19780000 & sample == "Rep01_W_F" & method == "adapt")

if (nrow(prod_row) == 0) {
  stop("Position 19780000, sample Rep01_W_F, method adapt not found in production data")
}

cat("Found production row for position 19780000, sample Rep01_W_F\n")

# Compare results
cat("\n=== COMPARISON ===\n")

# Compare haplotype estimates
hunk_haps <- hunk_result$Haps
prod_haps <- prod_row$Haps[[1]]

cat("Haplotype estimates:\n")
cat("Hunk:  ", paste(round(hunk_haps, 6), collapse = ", "), "\n")
cat("Prod:  ", paste(round(prod_haps, 6), collapse = ", "), "\n")
hap_diff <- sum(abs(hunk_haps - prod_haps))
cat("Sum of absolute differences:", round(hap_diff, 6), "\n")

# Compare error matrices
hunk_err <- hunk_result$Err
prod_err <- prod_row$Err[[1]]

cat("\nError matrix diagonal sums:\n")
hunk_diag_sum <- sum(diag(hunk_err))
prod_diag_sum <- sum(diag(prod_err))
cat("Hunk diagonal sum:", round(hunk_diag_sum, 6), "\n")
cat("Prod diagonal sum:", round(prod_diag_sum, 6), "\n")
cat("Difference (hunk - prod):", round(hunk_diag_sum - prod_diag_sum, 6), "\n")
cat("Ratio (hunk / prod):", round(hunk_diag_sum / prod_diag_sum, 6), "\n")

# Compare groups
hunk_groups <- hunk_result$Groups
prod_groups <- prod_row$Groups[[1]]

cat("\nGroups:\n")
cat("Hunk:  ", paste(hunk_groups, collapse = ", "), "\n")
cat("Prod:  ", paste(prod_groups, collapse = ", "), "\n")
cat("Groups match:", identical(hunk_groups, prod_groups), "\n")

cat("\n=== COMPARISON COMPLETE ===\n")
