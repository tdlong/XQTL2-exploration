#!/usr/bin/env Rscript

# Simple Haplotype Check
# Answers 3 questions:
# 1. Is the analysis done (haplotype file exists)?
# 2. What fraction of positions have founder NAs vs estimates (per sample)?
# 3. What fraction of positions are estimate_OK?

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("=== SIMPLE HAPLOTYPE CHECK ===\n\n")

# Parameters
chr <- "chr2R"
results_dir <- "process/JUICE/haplotype_results"

# Load parameter file to get founders
source("helpfiles/JUICE_haplotype_parameters.R")
cat("Checking for", length(founders), "founders:", paste(founders, collapse = ", "), "\n\n")

# Expected haplotype files
methods <- c(
  "fixed_window_20kb", "fixed_window_50kb", "fixed_window_100kb", "fixed_window_200kb", "fixed_window_500kb",
  "adaptive_window_h4", "adaptive_window_h6", "adaptive_window_h8", "adaptive_window_h10"
)

for (method in methods) {
  file_path <- file.path(results_dir, paste0(method, "_results_", chr, ".RDS"))
  
  cat("=== METHOD:", method, "===\n")
  
  # Question 1: Is the analysis done?
  if (file.exists(file_path)) {
    cat("✓ ANALYSIS DONE: File exists\n")
    
    # Load the data
    results <- readRDS(file_path)
    cat("  Data:", nrow(results), "rows\n")
    cat("  Columns:", paste(names(results), collapse = ", "), "\n")
    
    # Question 2: What fraction of positions have estimates vs NAs?
    # In binary distinguishability format: estimate_OK tells us if founders could be distinguished
    cat("✓ Using binary distinguishability format (estimate_OK column)\n")
    cat("  estimate_OK = 1: All founders distinguishable (haplotype estimatable)\n")
    cat("  estimate_OK = 0: Founders not distinguishable (haplotype NOT estimatable)\n")
    
    # Question 3: What fraction are estimate_OK?
    if ("estimate_OK" %in% names(results)) {
      if ("sample" %in% names(results)) {
        # Per-sample estimate_OK
        ok_summary <- results %>%
          group_by(sample) %>%
          summarise(
            total_positions = n(),
            estimate_OK_1 = sum(estimate_OK == 1, na.rm = TRUE),
            estimate_OK_0 = sum(estimate_OK == 0, na.rm = TRUE),
            estimate_OK_NA = sum(is.na(estimate_OK)),
            .groups = "drop"
          ) %>%
          mutate(fraction_OK = round(estimate_OK_1 / total_positions, 3))
        
        cat("  Per-sample estimate_OK:\n")
        for (i in 1:nrow(ok_summary)) {
          cat(sprintf("    %s: %.3f estimate_OK (1=%d, 0=%d, NA=%d)\n",
                     ok_summary$sample[i],
                     ok_summary$fraction_OK[i],
                     ok_summary$estimate_OK_1[i],
                     ok_summary$estimate_OK_0[i],
                     ok_summary$estimate_OK_NA[i]))
        }
      } else {
        # Single estimate_OK
        total_pos <- nrow(results)
        ok_1 <- sum(results$estimate_OK == 1, na.rm = TRUE)
        ok_0 <- sum(results$estimate_OK == 0, na.rm = TRUE)
        ok_na <- sum(is.na(results$estimate_OK))
        
        cat(sprintf("  estimate_OK: %.3f (1=%d, 0=%d, NA=%d)\n",
                   ok_1/total_pos, ok_1, ok_0, ok_na))
      }
    } else {
      cat("✗ No estimate_OK column found\n")
    }
    
  } else {
    cat("✗ ANALYSIS NOT DONE: File missing\n")
  }
  
  cat("\n")
}

cat("=== CHECK COMPLETE ===\n")
