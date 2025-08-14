#!/usr/bin/env Rscript

# Investigate Haplotype Calling Process
# Understand why duplicates are being produced

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"

cat("=== Investigating Haplotype Calling Process ===\n")
cat("Chromosome:", chr, "\n\n")

# Load all haplotype results
haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
if (file.exists(haplotype_file)) {
  haplotype_results <- read_rds(haplotype_file)
  cat("✓ Loaded haplotype results:", nrow(haplotype_results), "rows\n")
} else {
  cat("❌ Haplotype file not found:", haplotype_file, "\n")
  quit(status = 1)
}

# Check structure
cat("\n=== Data Structure ===\n")
cat("Columns:", paste(names(haplotype_results), collapse = ", "), "\n")
cat("Unique h_cutoff values:", paste(sort(unique(haplotype_results$h_cutoff)), collapse = ", "), "\n")
cat("Unique samples:", paste(sort(unique(haplotype_results$sample)), collapse = ", "), "\n")
cat("Unique founders:", paste(sort(unique(haplotype_results$founder)), collapse = ", "), "\n")

# Check for duplicates across all estimators
cat("\n=== Duplicate Analysis Across All Estimators ===\n")
duplicate_summary <- haplotype_results %>%
  group_by(h_cutoff, sample, pos, founder) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(h_cutoff) %>%
  summarise(
    total_positions = n(),
    duplicate_positions = sum(n > 1),
    duplicate_percent = round(duplicate_positions / total_positions * 100, 1),
    max_duplicates = max(n),
    .groups = "drop"
  )

print(duplicate_summary)

# Check specific samples
cat("\n=== Sample-by-Sample Duplicate Check ===\n")
sample_duplicates <- haplotype_results %>%
  group_by(h_cutoff, sample, pos, founder) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(h_cutoff, sample) %>%
  summarise(
    total_positions = n(),
    duplicate_positions = sum(n > 1),
    duplicate_percent = round(duplicate_positions / total_positions * 100, 1),
    max_duplicates = max(n),
    .groups = "drop"
  ) %>%
  arrange(h_cutoff, desc(duplicate_percent))

print(sample_duplicates)

# Check specific positions with duplicates
cat("\n=== Examining Duplicate Positions ===\n")
duplicate_positions <- haplotype_results %>%
  group_by(h_cutoff, sample, pos, founder) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1) %>%
  arrange(h_cutoff, sample, pos, founder)

if (nrow(duplicate_positions) > 0) {
  cat("First 10 duplicate positions:\n")
  print(head(duplicate_positions, 10))
  
  # Look at a specific duplicate case
  cat("\n=== Detailed Look at First Duplicate ===\n")
  first_dup <- duplicate_positions[1, ]
  detailed_dup <- haplotype_results %>%
    filter(h_cutoff == first_dup$h_cutoff,
           sample == first_dup$sample,
           pos == first_dup$pos,
           founder == first_dup$founder)
  
  cat("Sample:", first_dup$sample, "Position:", first_dup$pos, "Founder:", first_dup$founder, "\n")
  cat("Number of duplicates:", nrow(detailed_dup), "\n")
  cat("All frequencies:", paste(detailed_dup$freq, collapse = ", "), "\n")
  
  # Check if there are other columns that differ
  cat("\nChecking for differences in other columns:\n")
  for (col in names(detailed_dup)) {
    if (col != "freq") {
      unique_vals <- unique(detailed_dup[[col]])
      if (length(unique_vals) > 1) {
        cat("  ", col, ":", paste(unique_vals, collapse = ", "), "\n")
      }
    }
  }
}

# Check if this is happening in other chromosomes
cat("\n=== Checking Other Chromosomes ===\n")
other_chr_files <- list.files(output_dir, pattern = "adaptive_window_results_.*\\.RDS", full.names = TRUE)
other_chr_files <- other_chr_files[!grepl(chr, other_chr_files)]

if (length(other_chr_files) > 0) {
  cat("Found other chromosome files:", length(other_chr_files), "\n")
  
  # Check first other chromosome
  other_chr <- other_chr_files[1]
  cat("Checking:", basename(other_chr), "\n")
  
  other_results <- read_rds(other_chr)
  other_duplicates <- other_results %>%
    group_by(h_cutoff, sample, pos, founder) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)
  
  cat("Duplicates in other chromosome:", nrow(other_duplicates), "\n")
  if (nrow(other_duplicates) > 0) {
    cat("Max duplicates per founder:", max(other_duplicates$n), "\n")
  }
} else {
  cat("No other chromosome files found\n")
}

cat("\n=== Investigation Complete ===\n")
cat("Next steps:\n")
cat("1. Check if this is happening in fixed vs adaptive estimators\n")
cat("2. Look at the haplotype calling code to understand why duplicates are created\n")
cat("3. Check if this is a data processing issue or algorithm issue\n")
