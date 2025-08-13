#!/usr/bin/env Rscript

# Quick SNP Coverage Check
# Count how many SNPs have valid intervals

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"
sample_name <- "GJ_3_1"

cat("=== SNP Coverage Check ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n\n")

# Load haplotype results
h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
haplotype_results <- read_rds(haplotype_file) %>%
  filter(h_cutoff == !!h_cutoff)

# Filter to sample and euchromatin
euchromatin_start <- 5398184
euchromatin_end <- 24684540

sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name, 
         pos >= euchromatin_start, 
         pos <= euchromatin_end)

cat("✓ Sample haplotypes loaded:", nrow(sample_haplotypes), "rows\n")
cat("Position range:", min(sample_haplotypes$pos), "-", max(sample_haplotypes$pos), "bp\n\n")

# Load SNP data
refalt_file <- file.path(output_dir, paste0("df3.", chr, ".RDS"))
df2 <- read_rds(refalt_file)

# Filter SNPs to euchromatin
good_snps <- df2 %>%
  filter(name %in% c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end)

cat("✓ Valid euchromatic SNPs:", nrow(valid_snps), "\n")
cat("SNP position range:", min(valid_snps$POS), "-", max(valid_snps$POS), "bp\n\n")

# Get haplotype positions
haplotype_positions <- sort(unique(sample_haplotypes$pos))

# Check coverage for a sample of SNPs
cat("=== Checking SNP Coverage ===\n")
cat("Testing first 1000 SNPs...\n")

test_snps <- head(valid_snps$POS, 1000)
valid_intervals <- 0
edge_cases <- 0
no_coverage <- 0

for (snp_pos in test_snps) {
  # Find haplotype positions <= snp_pos (left side)
  left_positions <- haplotype_positions[haplotype_positions <= snp_pos]
  left_pos <- if (length(left_positions) > 0) max(left_positions) else NA
  
  # Find haplotype positions > snp_pos (right side)  
  right_positions <- haplotype_positions[haplotype_positions > snp_pos]
  right_pos <- if (length(right_positions) > 0) min(right_positions) else NA
  
  if (!is.na(left_pos) && !is.na(right_pos)) {
    valid_intervals <- valid_intervals + 1
  } else if (is.na(left_pos) && is.na(right_pos)) {
    no_coverage <- no_coverage + 1
  } else {
    edge_cases <- edge_cases + 1
  }
}

cat("Results for first 1000 SNPs:\n")
cat("  Valid intervals (both sides):", valid_intervals, "(", round(valid_intervals/1000*100, 1), "%)\n")
cat("  Edge cases (one side only):", edge_cases, "(", round(edge_cases/1000*100, 1), "%)\n")
cat("  No coverage:", no_coverage, "(", round(no_coverage/1000*100, 1), "%)\n\n")

# Check middle section
cat("Testing middle 1000 SNPs...\n")
mid_start <- length(valid_snps$POS) %/% 2 - 500
mid_snps <- valid_snps$POS[mid_start:(mid_start + 999)]

valid_intervals_mid <- 0
edge_cases_mid <- 0
no_coverage_mid <- 0

for (snp_pos in mid_snps) {
  left_positions <- haplotype_positions[haplotype_positions <= snp_pos]
  left_pos <- if (length(left_positions) > 0) max(left_positions) else NA
  
  right_positions <- haplotype_positions[haplotype_positions > snp_pos]
  right_pos <- if (length(right_positions) > 0) min(right_positions) else NA
  
  if (!is.na(left_pos) && !is.na(right_pos)) {
    valid_intervals_mid <- valid_intervals_mid + 1
  } else if (is.na(left_pos) && is.na(right_pos)) {
    no_coverage_mid <- no_coverage_mid + 1
  } else {
    edge_cases_mid <- edge_cases_mid + 1
  }
}

cat("Results for middle 1000 SNPs:\n")
cat("  Valid intervals (both sides):", valid_intervals_mid, "(", round(valid_intervals_mid/1000*100, 1), "%)\n")
cat("  Edge cases (one side only):", edge_cases_mid, "(", round(edge_cases_mid/1000*100, 1), "%)\n")
cat("  No coverage:", no_coverage_mid, "(", round(no_coverage_mid/1000*100, 1), "%)\n\n")

cat("=== Coverage Check Complete ===\n")
