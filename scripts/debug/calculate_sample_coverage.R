#!/usr/bin/env Rscript

# Calculate average coverage for specific samples from RefAlt files
# This script analyzes REF + ALT allele counts to determine sequencing coverage

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript calculate_sample_coverage.R <chr> <juice_param_file> <zinc2_param_file>")
}

chr <- args[1]
juice_param_file <- args[2]
zinc2_param_file <- args[3]

cat("=== CALCULATING SAMPLE COVERAGE ===\n")
cat("Chromosome:", chr, "\n")
cat("JUICE parameter file:", juice_param_file, "\n")
cat("ZINC2 parameter file:", zinc2_param_file, "\n\n")

# Load parameters
source(juice_param_file)
juice_refalt_file <- file.path("process/JUICE", paste0("RefAlt.", chr, ".txt"))

source(zinc2_param_file)
zinc2_refalt_file <- file.path("process/ZINC2", paste0("RefAlt.", chr, ".txt"))

# Function to calculate coverage for a sample
calculate_sample_coverage <- function(refalt_file, sample_name, dataset_name) {
  cat("=== ANALYZING", toupper(dataset_name), "===\n")
  
  if (!file.exists(refalt_file)) {
    cat("❌ RefAlt file not found:", refalt_file, "\n")
    return(NULL)
  }
  
  cat("✓ Loading RefAlt file...\n")
  refalt_data <- read.table(refalt_file, header = TRUE)
  
  cat("RefAlt data dimensions:", nrow(refalt_data), "rows,", ncol(refalt_data), "columns\n")
  cat("Columns:", paste(names(refalt_data), collapse = ", "), "\n\n")
  
  # Find columns for the specific sample
  ref_col <- paste0("REF_", sample_name)
  alt_col <- paste0("ALT_", sample_name)
  
  if (!ref_col %in% names(refalt_data)) {
    cat("❌ REF column not found:", ref_col, "\n")
    cat("Available REF columns:", paste(grep("^REF_", names(refalt_data), value = TRUE), collapse = ", "), "\n")
    return(NULL)
  }
  
  if (!alt_col %in% names(refalt_data)) {
    cat("❌ ALT column not found:", alt_col, "\n")
    cat("Available ALT columns:", paste(grep("^ALT_", names(refalt_data), value = TRUE), collapse = ", "), "\n")
    return(NULL)
  }
  
  cat("✓ Found columns:", ref_col, "and", alt_col, "\n")
  
  # Calculate coverage (REF + ALT)
  coverage_data <- refalt_data %>%
    select(CHROM, POS, all_of(c(ref_col, alt_col))) %>%
    mutate(
      ref_count = !!sym(ref_col),
      alt_count = !!sym(alt_col),
      coverage = ref_count + alt_count
    ) %>%
    filter(!is.na(coverage))  # Remove any rows with missing data
  
  cat("SNPs with coverage data:", nrow(coverage_data), "\n")
  
  # Calculate coverage statistics
  coverage_stats <- coverage_data %>%
    summarize(
      total_snps = n(),
      mean_coverage = mean(coverage, na.rm = TRUE),
      median_coverage = median(coverage, na.rm = TRUE),
      min_coverage = min(coverage, na.rm = TRUE),
      max_coverage = max(coverage, na.rm = TRUE),
      sd_coverage = sd(coverage, na.rm = TRUE),
      snps_with_coverage_0 = sum(coverage == 0),
      snps_with_coverage_1_10 = sum(coverage >= 1 & coverage <= 10),
      snps_with_coverage_11_50 = sum(coverage >= 11 & coverage <= 50),
      snps_with_coverage_51_100 = sum(coverage >= 51 & coverage <= 100),
      snps_with_coverage_over_100 = sum(coverage > 100),
      .groups = "drop"
    )
  
  cat("\nCoverage Statistics for", sample_name, ":\n")
  cat("  Total SNPs:", coverage_stats$total_snps, "\n")
  cat("  Mean coverage:", round(coverage_stats$mean_coverage, 2), "\n")
  cat("  Median coverage:", coverage_stats$median_coverage, "\n")
  cat("  Coverage range:", coverage_stats$min_coverage, "-", coverage_stats$max_coverage, "\n")
  cat("  Standard deviation:", round(coverage_stats$sd_coverage, 2), "\n")
  
  cat("\nCoverage Distribution:\n")
  cat("  Coverage = 0:", coverage_stats$snps_with_coverage_0, "SNPs\n")
  cat("  Coverage 1-10:", coverage_stats$snps_with_coverage_1_10, "SNPs\n")
  cat("  Coverage 11-50:", coverage_stats$snps_with_coverage_11_50, "SNPs\n")
  cat("  Coverage 51-100:", coverage_stats$snps_with_coverage_51_100, "SNPs\n")
  cat("  Coverage > 100:", coverage_stats$snps_with_coverage_over_100, "SNPs\n")
  
  # Calculate coverage percentiles
  coverage_percentiles <- quantile(coverage_data$coverage, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)
  cat("\nCoverage Percentiles:\n")
  for (i in 1:length(coverage_percentiles)) {
    pct_name <- names(coverage_percentiles)[i]
    pct_value <- coverage_percentiles[i]
    cat(sprintf("  %s: %.1f\n", pct_name, pct_value))
  }
  
  cat("\n")
  
  return(coverage_stats)
}

# Calculate coverage for JUICE sample
juice_coverage <- calculate_sample_coverage(juice_refalt_file, "AJ_1_1", "JUICE")

# Calculate coverage for ZINC2 sample
zinc2_coverage <- calculate_sample_coverage(zinc2_refalt_file, "Rep01_W_F", "ZINC2")

# Compare the two samples
cat("=== COVERAGE COMPARISON ===\n")

if (!is.null(juice_coverage) && !is.null(zinc2_coverage)) {
  cat("Mean Coverage Comparison:\n")
  cat("  JUICE (AJ_1_1):", round(juice_coverage$mean_coverage, 2), "\n")
  cat("  ZINC2 (Rep01_W_F):", round(zinc2_coverage$mean_coverage, 2), "\n")
  
  coverage_ratio <- juice_coverage$mean_coverage / zinc2_coverage$mean_coverage
  cat("  Coverage ratio (JUICE/ZINC2):", round(coverage_ratio, 2), "\n")
  
  if (coverage_ratio > 1.1) {
    cat("  → JUICE has", round((coverage_ratio - 1) * 100, 1), "% higher coverage\n")
  } else if (coverage_ratio < 0.9) {
    cat("  → ZINC2 has", round((1/coverage_ratio - 1) * 100, 1), "% higher coverage\n")
  } else {
    cat("  → Similar coverage levels\n")
  }
  
  cat("\nMedian Coverage Comparison:\n")
  cat("  JUICE (AJ_1_1):", juice_coverage$median_coverage, "\n")
  cat("  ZINC2 (Rep01_W_F):", zinc2_coverage$median_coverage, "\n")
  
  cat("\nTotal SNPs Comparison:\n")
  cat("  JUICE (AJ_1_1):", juice_coverage$total_snps, "SNPs\n")
  cat("  ZINC2 (Rep01_W_F):", zinc2_coverage$total_snps, "SNPs\n")
  
} else {
  cat("❌ Could not compare coverage - one or both samples failed to load\n")
}

cat("\n=== COVERAGE ANALYSIS COMPLETE ===\n")
cat("This shows the sequencing depth differences between JUICE and ZINC2 samples.\n")
