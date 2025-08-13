#!/usr/bin/env Rscript

# Debug Streaming Algorithm Script
# Tests the interval finding logic step-by-step

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"
sample_name <- "GJ_3_1"

cat("=== Debug Streaming Algorithm ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n\n")

# Load haplotype results
cat("Loading haplotype results...\n")
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
cat("Loading SNP data...\n")
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
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(df2, multiple = "all")

cat("✓ Valid euchromatic SNPs:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n")
cat("SNP position range:", min(valid_snps$POS), "-", max(valid_snps$POS), "bp\n\n")

# Test the streaming algorithm logic
cat("=== Testing Streaming Algorithm Logic ===\n")

# Convert haplotypes to wide format
haplotype_freqs <- sample_haplotypes %>%
  pivot_wider(names_from = founder, values_from = freq, values_fill = NA)

# Get unique haplotype positions (sorted)
haplotype_positions <- sort(unique(haplotype_freqs$pos))
cat("Haplotype positions:", length(haplotype_positions), "\n")
cat("First 5:", paste(head(haplotype_positions, 5), collapse = ", "), "\n")
cat("Last 5:", paste(tail(haplotype_positions, 5), collapse = ", "), "\n\n")

# Get SNP positions to test
snp_positions <- valid_snps %>% distinct(POS) %>% pull(POS) %>% sort()
cat("SNP positions:", length(snp_positions), "\n")
cat("First 5:", paste(head(snp_positions, 5), collapse = ", "), "\n")
cat("Last 5:", paste(tail(snp_positions, 5), collapse = ", "), "\n\n")

# Test interval finding for first 10 SNPs
cat("=== Testing Interval Finding for First 10 SNPs ===\n")
test_snps <- head(snp_positions, 10)

for (i in seq_along(test_snps)) {
  snp_pos <- test_snps[i]
  cat("SNP", i, "at position", snp_pos, ":\n")
  
  # Find haplotype positions <= snp_pos (left side)
  left_positions <- haplotype_positions[haplotype_positions <= snp_pos]
  left_pos <- if (length(left_positions) > 0) max(left_positions) else NA
  
  # Find haplotype positions > snp_pos (right side)  
  right_positions <- haplotype_positions[haplotype_positions > snp_pos]
  right_pos <- if (length(right_positions) > 0) min(right_positions) else NA
  
  cat("  Left haplotype:", ifelse(is.na(left_pos), "NA", left_pos), "\n")
  cat("  Right haplotype:", ifelse(is.na(right_pos), "NA", right_pos), "\n")
  
  if (!is.na(left_pos) && !is.na(right_pos)) {
    cat("  ✓ Valid interval found\n")
  } else if (is.na(left_pos) && is.na(right_pos)) {
    cat("  ❌ No valid interval\n")
  } else {
    cat("  ⚠️  Edge case (one side only)\n")
  }
  cat("\n")
}

# Test the streaming algorithm's interval tracking
cat("=== Testing Streaming Algorithm Interval Tracking ===\n")
current_left_idx <- 1
current_right_idx <- 2

cat("Starting indices: left =", current_left_idx, ", right =", current_right_idx, "\n")

for (i in 1:5) {
  snp_pos <- test_snps[i]
  cat("SNP", i, "at position", snp_pos, ":\n")
  
  # Find the haplotype interval containing this SNP
  while (current_right_idx <= length(haplotype_positions) && 
         haplotype_positions[current_right_idx] < snp_pos) {
    current_left_idx <- current_right_idx
    current_right_idx <- current_right_idx + 1
    cat("  Advanced interval: left =", current_left_idx, ", right =", current_right_idx, "\n")
  }
  
  # Check if we have a valid interval
  if (current_left_idx > length(haplotype_positions) || 
      current_right_idx > length(haplotype_positions)) {
    cat("  ❌ Invalid interval indices\n")
    break
  }
  
  left_pos <- haplotype_positions[current_left_idx]
  right_pos <- haplotype_positions[current_right_idx]
  
  cat("  Final interval: [", left_pos, ",", right_pos, "]\n")
  cat("  SNP position:", snp_pos, "is between haplotypes\n")
  cat("\n")
}

cat("=== Debug Complete ===\n")
