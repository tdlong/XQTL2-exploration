#!/usr/bin/env Rscript

# Check Haplotype Flag
# Simple script to check what flag was set for a specific position in haplotype results
#
# Usage: Rscript check_haplotype_flag.R <chr> <output_dir> <estimator> <position>

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript check_haplotype_flag.R <chr> <output_dir> <estimator> <position>")
}

chr <- args[1]
output_dir <- args[2]
estimator <- args[3]
target_position <- as.numeric(args[4])

cat("=== Check Haplotype Flag ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n")
cat("Estimator:", estimator, "\n")
cat("Target position:", target_position, "\n\n")

# Load haplotype results
cat("Loading haplotype results...\n")
if (grepl("^fixed_", estimator)) {
  window_size_kb <- as.numeric(gsub("fixed_", "", gsub("kb", "", estimator)))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
} else if (grepl("^adaptive_h", estimator)) {
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
} else {
  stop("Invalid estimator format")
}

if (!file.exists(haplotype_file)) {
  stop("Haplotype file not found:", haplotype_file)
}

haplotype_results <- read_rds(haplotype_file)
cat("✓ Haplotype results loaded:", nrow(haplotype_results), "rows\n")

# Get founders from haplotype results
founders <- names(haplotype_results)[grepl("^[A-Z][0-9]+$|^AB[0-9]+$", names(haplotype_results))]
cat("Founders detected:", paste(founders, collapse = ", "), "\n")

# Check for the specific position
cat("\n=== Checking Position", target_position, "===\n")
position_results <- haplotype_results %>%
  filter(pos == target_position)

if (nrow(position_results) == 0) {
  cat("❌ No results found for position", target_position, "\n")
} else {
  cat("✓ Found", nrow(position_results), "results for position", target_position, "\n")
  
  # Show results for each sample
  for (sample_name in unique(position_results$sample)) {
    cat("\n--- Sample:", sample_name, "---\n")
    
    sample_result <- position_results %>%
      filter(sample == sample_name)
    
    if (nrow(sample_result) == 0) {
      cat("  No data for this sample\n")
    } else {
      # Check if estimate_OK column exists
      if ("estimate_OK" %in% names(sample_result)) {
        estimate_ok <- sample_result %>% pull(estimate_OK)
        cat("  estimate_OK flag:", estimate_ok, "\n")
        cat("  Flag interpretation:", ifelse(estimate_ok == 1, "RELIABLE", ifelse(estimate_ok == 0, "UNRELIABLE", "FAILED")), "\n")
      } else {
        cat("  No estimate_OK column found\n")
      }
      
      # Show haplotype frequencies
      founder_cols <- intersect(founders, names(sample_result))
      if (length(founder_cols) > 0) {
        cat("  Haplotype frequencies:\n")
        for (founder in founder_cols) {
          freq <- sample_result %>% pull(founder)
          cat("    ", founder, ":", freq, "\n")
        }
      } else {
        cat("  No founder frequency columns found\n")
      }
      
      # Show all available columns
      cat("  Available columns:", paste(names(sample_result), collapse = ", "), "\n")
    }
  }
}

cat("\n=== Check Complete ===\n")
cat("This shows what flag was actually set by the production algorithm.\n")
cat("Check if estimate_OK = 0 (unreliable) when founders cannot be distinguished.\n")
