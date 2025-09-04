#!/usr/bin/env Rscript

# Compare region 870 between JUICE and ZINC2 to identify the MAE plotting issue
# This focuses specifically on the region being plotted

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0) {
  stop("Usage: Rscript compare_region_870_juice_vs_zinc2.R")
}

cat("=== COMPARING REGION 870: JUICE vs ZINC2 ===\n")
cat("Region: 8,700,000 ± 100,000 bp (position 870 in 10kb units)\n\n")

# Define the region
midpoint_10kb <- 870
region_start_bp <- midpoint_10kb * 10000 - 100000  # 100kb on each side
region_end_bp <- midpoint_10kb * 10000 + 100000

cat("Region boundaries:", region_start_bp, "-", region_end_bp, "bp\n")
cat("This corresponds to:", round(region_start_bp/10000), "-", round(region_end_bp/10000), "(10kb units)\n\n")

# Function to analyze a dataset
analyze_dataset <- function(dataset_name, param_file, output_dir) {
  cat("=== ANALYZING", toupper(dataset_name), "===\n")
  
  # Load parameters
  source(param_file)
  results_dir <- file.path(output_dir, "haplotype_results")
  
  # Load summary file
  summary_file <- file.path(results_dir, paste0("summary_chr2R.RDS"))
  if (!file.exists(summary_file)) {
    cat("❌ Summary file not found:", summary_file, "\n")
    return(NULL)
  }
  
  cat("✓ Loading summary file...\n")
  summary_data <- readRDS(summary_file)
  
  # Get first sample
  first_sample <- unique(summary_data$sample)[1]
  cat("First sample:", first_sample, "\n")
  
  # Filter to region and first sample
  region_data <- summary_data %>%
    filter(sample == first_sample) %>%
    filter(pos >= region_start_bp & pos <= region_end_bp) %>%
    filter(method %in% c("adaptive_h4", "fixed_20kb", "fixed_100kb", "smooth_h4"))
  
  cat("Data points in region:", nrow(region_data), "\n")
  cat("Methods found:", paste(unique(region_data$method), collapse = ", "), "\n\n")
  
  if (nrow(region_data) == 0) {
    cat("❌ No data found in region!\n")
    return(NULL)
  }
  
  # Analyze MAE data by method
  cat("MAE Analysis by Method:\n")
  mae_analysis <- region_data %>%
    group_by(method) %>%
    summarize(
      total_positions = n(),
      positions_with_mae = sum(!is.na(MAE)),
      positions_without_mae = sum(is.na(MAE)),
      mae_coverage_pct = round(positions_with_mae / total_positions * 100, 1),
      mean_mae = mean(MAE, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(mae_analysis)
  cat("\n")
  
  # Show detailed data for each method
  for (method in unique(region_data$method)) {
    method_data <- region_data %>% filter(method == method)
    cat("=== METHOD:", method, "===\n")
    cat("Total positions:", nrow(method_data), "\n")
    cat("Positions with MAE:", sum(!is.na(method_data$MAE)), "\n")
    cat("Positions without MAE:", sum(is.na(method_data$MAE)), "\n")
    
    if (sum(!is.na(method_data$MAE)) > 0) {
      cat("MAE range:", round(range(method_data$MAE, na.rm = TRUE), 4), "\n")
      cat("MAE mean:", round(mean(method_data$MAE, na.rm = TRUE), 4), "\n")
    }
    
    # Show first few rows
    cat("First 5 positions:\n")
    print(method_data %>% head(5) %>% select(pos, MAE, estimate_OK, NSNPs))
    cat("\n")
  }
  
  return(region_data)
}

# Analyze JUICE
juice_data <- analyze_dataset("JUICE", "helpfiles/JUICE_haplotype_parameters.R", "process/JUICE")

# Analyze ZINC2
zinc2_data <- analyze_dataset("ZINC2", "helpfiles/ZINC2_haplotype_parameters.R", "process/ZINC2")

# Compare the datasets
cat("=== COMPARISON SUMMARY ===\n")

if (!is.null(juice_data) && !is.null(zinc2_data)) {
  cat("JUICE data points:", nrow(juice_data), "\n")
  cat("ZINC2 data points:", nrow(zinc2_data), "\n")
  
  # Compare MAE coverage
  juice_mae <- juice_data %>%
    group_by(method) %>%
    summarize(juice_mae_coverage = sum(!is.na(MAE)) / n(), .groups = "drop")
  
  zinc2_mae <- zinc2_data %>%
    group_by(method) %>%
    summarize(zinc2_mae_coverage = sum(!is.na(MAE)) / n(), .groups = "drop")
  
  comparison <- juice_mae %>%
    left_join(zinc2_mae, by = "method") %>%
    mutate(
      juice_mae_coverage = round(juice_mae_coverage * 100, 1),
      zinc2_mae_coverage = round(zinc2_mae_coverage * 100, 1)
    )
  
  cat("\nMAE Coverage Comparison:\n")
  print(comparison)
  
  # Check position ranges
  cat("\nPosition Ranges:\n")
  cat("JUICE positions:", range(juice_data$pos), "\n")
  cat("ZINC2 positions:", range(zinc2_data$pos), "\n")
  
} else {
  cat("❌ Could not compare datasets - one or both failed to load\n")
}

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("This should reveal why JUICE shows missing MAE data in the region plot.\n")
