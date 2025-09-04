#!/usr/bin/env Rscript

# Diagnostic script to identify why MAE data is missing for JUICE
# Run this on the cluster to debug the MAE calculation issue

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript diagnose_juice_mae_issue.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

cat("=== DIAGNOSING JUICE MAE ISSUE ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

# Load the summary file to see what's there
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
if (!file.exists(summary_file)) {
  stop("Summary file not found: ", summary_file)
}

cat("Loading summary file...\n")
summary_data <- readRDS(summary_file)

cat("✓ Summary data loaded:", nrow(summary_data), "rows\n")
cat("Methods available:", paste(unique(summary_data$method), collapse = ", "), "\n")
cat("Samples available:", paste(unique(summary_data$sample), collapse = ", "), "\n\n")

# Check MAE data availability by method
cat("=== MAE DATA AVAILABILITY BY METHOD ===\n")
mae_summary <- summary_data %>%
  group_by(method) %>%
  summarize(
    total_positions = n(),
    positions_with_mae = sum(!is.na(MAE)),
    positions_without_mae = sum(is.na(MAE)),
    mae_coverage_pct = round(positions_with_mae / total_positions * 100, 1),
    .groups = "drop"
  )

print(mae_summary)
cat("\n")

# Focus on adaptive_h4 for detailed analysis
cat("=== DETAILED ANALYSIS FOR adaptive_h4 ===\n")
adaptive_data <- summary_data %>%
  filter(method == "adaptive_h4")

if (nrow(adaptive_data) > 0) {
  first_sample <- unique(adaptive_data$sample)[1]
  cat("Using first sample:", first_sample, "\n\n")
  
  sample_data <- adaptive_data %>%
    filter(sample == first_sample)
  
  cat("Sample data points:", nrow(sample_data), "\n")
  cat("Positions with MAE:", sum(!is.na(sample_data$MAE)), "\n")
  cat("Positions without MAE:", sum(is.na(sample_data$MAE)), "\n")
  
  # Check the actual MAE values
  mae_values <- sample_data$MAE[!is.na(sample_data$MAE)]
  if (length(mae_values) > 0) {
    cat("MAE value range:", round(range(mae_values), 4), "\n")
    cat("MAE mean:", round(mean(mae_values), 4), "\n")
  } else {
    cat("No MAE values found!\n")
  }
  
  # Show first few rows with and without MAE
  cat("\nFirst 5 positions WITH MAE:\n")
  print(sample_data %>% filter(!is.na(MAE)) %>% head(5) %>% select(pos, MAE, estimate_OK, NSNPs))
  
  cat("\nFirst 5 positions WITHOUT MAE:\n")
  print(sample_data %>% filter(is.na(MAE)) %>% head(5) %>% select(pos, MAE, estimate_OK, NSNPs))
  
} else {
  cat("No adaptive_h4 data found!\n")
}

# Now let's check the raw SNP imputation data
cat("\n=== CHECKING RAW SNP IMPUTATION DATA ===\n")
snp_file <- file.path(results_dir, paste0("snp_imputation_adaptive_h4_", chr, ".RDS"))
if (!file.exists(snp_file)) {
  cat("❌ SNP imputation file not found:", snp_file, "\n")
} else {
  cat("✓ Loading SNP imputation data...\n")
  snp_data <- readRDS(snp_file)
  
  cat("SNP data dimensions:", nrow(snp_data), "rows,", ncol(snp_data), "columns\n")
  cat("Columns:", paste(names(snp_data), collapse = ", "), "\n")
  cat("Samples in SNP data:", paste(unique(snp_data$sample), collapse = ", "), "\n")
  
  # Check for NA values in imputed column
  cat("\nImputed column analysis:\n")
  cat("Total SNPs:", nrow(snp_data), "\n")
  cat("SNPs with imputed = NA:", sum(is.na(snp_data$imputed)), "\n")
  cat("SNPs with imputed != NA:", sum(!is.na(snp_data$imputed)), "\n")
  
  # Check position range
  cat("\nPosition range in SNP data:", range(snp_data$pos), "\n")
  
  # Check a specific sample
  if (nrow(snp_data) > 0) {
    first_sample_snp <- unique(snp_data$sample)[1]
    cat("First sample in SNP data:", first_sample_snp, "\n")
    
    sample_snp_data <- snp_data %>%
      filter(sample == first_sample_snp)
    
    cat("SNPs for first sample:", nrow(sample_snp_data), "\n")
    cat("SNPs with imputed != NA for first sample:", sum(!is.na(sample_snp_data$imputed)), "\n")
    
    # Show first few rows
    cat("\nFirst 5 rows of SNP data for first sample:\n")
    print(sample_snp_data %>% head(5) %>% select(pos, observed, imputed, sample))
  }
}

# Check haplotype data
cat("\n=== CHECKING HAPLOTYPE DATA ===\n")
haplo_file <- file.path(results_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
if (!file.exists(haplo_file)) {
  cat("❌ Haplotype file not found:", haplo_file, "\n")
} else {
  cat("✓ Loading haplotype data...\n")
  haplo_data <- readRDS(haplo_file)
  
  cat("Haplotype data dimensions:", nrow(haplo_data), "rows,", ncol(haplo_data), "columns\n")
  cat("Columns:", paste(names(haplo_data), collapse = ", "), "\n")
  cat("Samples in haplotype data:", paste(unique(haplo_data$sample), collapse = ", "), "\n")
  
  # Check position range
  cat("\nPosition range in haplotype data:", range(haplo_data$pos), "\n")
  
  # Check estimate_OK distribution
  cat("\nEstimate_OK distribution:\n")
  print(table(haplo_data$estimate_OK, useNA = "ifany"))
}

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("This should help identify why MAE data is missing for JUICE.\n")
cat("Look for:\n")
cat("1. Sample name mismatches between SNP and haplotype data\n")
cat("2. All imputed values being NA in SNP data\n")
cat("3. No SNPs within ±5kb of haplotype positions\n")
cat("4. Position range mismatches between datasets\n")
