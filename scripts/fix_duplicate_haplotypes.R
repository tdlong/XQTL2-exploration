#!/usr/bin/env Rscript

# Fix Duplicate Haplotypes
# Remove duplicate founder frequencies and keep one per position

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"

cat("=== Fixing Duplicate Haplotypes ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n\n")

# Load haplotype results
h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
haplotype_results <- read_rds(haplotype_file) %>%
  filter(h_cutoff == !!h_cutoff)

cat("Original haplotype results:", nrow(haplotype_results), "rows\n")

# Check for duplicates
duplicate_check <- haplotype_results %>%
  group_by(sample, pos, founder) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

cat("Positions with duplicate founders:", nrow(duplicate_check), "\n")
if (nrow(duplicate_check) > 0) {
  cat("Sample with most duplicates:", duplicate_check$sample[which.max(duplicate_check$n)], "\n")
  cat("Max duplicates per founder:", max(duplicate_check$n), "\n\n")
}

# Fix duplicates by keeping the first occurrence of each founder per position
cat("Fixing duplicates...\n")
haplotype_fixed <- haplotype_results %>%
  group_by(sample, pos, founder) %>%
  slice_head(n = 1) %>%
  ungroup()

cat("Fixed haplotype results:", nrow(haplotype_fixed), "rows\n")
cat("Removed", nrow(haplotype_results) - nrow(haplotype_fixed), "duplicate rows\n\n")

# Verify fix
duplicate_check_after <- haplotype_fixed %>%
  group_by(sample, pos, founder) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if (nrow(duplicate_check_after) == 0) {
  cat("✓ All duplicates fixed!\n")
} else {
  cat("❌ Still have duplicates:", nrow(duplicate_check_after), "\n")
}

# Check specific positions
test_positions <- c(5400000, 5410000)
cat("\n=== Checking Fixed Positions ===\n")

for (pos in test_positions) {
  pos_data <- haplotype_fixed %>% filter(pos == !!pos, sample == "GJ_3_1")
  cat("Position", pos, "has", nrow(pos_data), "rows\n")
  
  if (nrow(pos_data) > 0) {
    cat("Founders:", paste(pos_data$founder, collapse = ", "), "\n")
    cat("Frequencies:", paste(round(pos_data$freq, 6), collapse = ", "), "\n")
  }
  cat("\n")
}

# Save fixed data
output_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, "_fixed.RDS"))
cat("Saving fixed data to:", output_file, "\n")
write_rds(haplotype_fixed, output_file)

cat("=== Fix Complete ===\n")
cat("Original file:", haplotype_file, "\n")
cat("Fixed file:", output_file, "\n")
cat("You can now use the fixed file for interpolation.\n")
