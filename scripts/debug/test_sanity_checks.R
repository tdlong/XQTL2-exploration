#!/usr/bin/env Rscript

# =============================================================================
# Test Script for Sanity Checks
# =============================================================================
# Tests the exact case that was failing to verify our sanity checks work
# Position: 510000, Sample: GJ_3_1, Window: 20kb
# Usage: Rscript scripts/test_sanity_checks.R

library(tidyverse)

# Set up the test case
test_chr <- "chr2R"
test_pos <- 510000
test_sample <- "GJ_3_1"
test_window_size <- 20000
test_mydir <- "process/JUICE"
test_parfile <- "helpfiles/JUICE/JUICE_haplotype_parameters.R"

cat("=== Testing Sanity Checks ===\n")
cat("Chromosome:", test_chr, "\n")
cat("Position:", test_pos, "\n")
cat("Sample:", test_sample, "\n")
cat("Window size:", test_window_size, "bp\n")
cat("Parameter file:", test_parfile, "\n")
cat("Data directory:", test_mydir, "\n\n")

# Source the parameter file
if (file.exists(test_parfile)) {
  source(test_parfile)
  cat("✓ Parameter file loaded\n")
  cat("Founders:", paste(founders, collapse = ", "), "\n")
  cat("h_cutoff:", h_cutoff, "\n\n")
} else {
  stop("Parameter file not found: ", test_parfile)
}

# Load the data
cat("Loading data...\n")
refalt_file <- paste0(test_mydir, "/RefAlt.", test_chr, ".txt")
if (!file.exists(refalt_file)) {
  stop("REFALT file not found: ", refalt_file)
}

df <- read.table(refalt_file, header = TRUE)
cat("✓ REFALT data loaded:", nrow(df), "rows,", ncol(df), "columns\n")

# Transform REF/ALT counts to frequencies
cat("Transforming data...\n")
df2 <- df %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    N = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

cat("✓ Data transformed\n")

# Filter for high-quality SNPs
cat("Filtering for high-quality SNPs...\n")
good_snps <- df2 %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

cat("✓ High-quality SNPs identified:", nrow(good_snps), "\n\n")

# Subset dataset to only high-quality SNPs
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Create window around test position
window_start <- max(0, test_pos - test_window_size/2)
window_end <- test_pos + test_window_size/2

cat("=== Window Analysis ===\n")
cat("Window:", window_start, "-", window_end, "bp\n")
cat("Window size:", test_window_size, "bp\n")

# Filter SNPs in this window
window_snps <- df3 %>%
  filter(CHROM == test_chr &
         POS > window_start &
         POS < window_end &
         (name %in% founders | name == test_sample)) %>%
  select(-c(CHROM, N)) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("✓ Window SNPs processed\n")
cat("SNPs in window:", nrow(window_snps), "\n")

# Check if we have the test sample
if (!(test_sample %in% names(window_snps))) {
  stop("Test sample", test_sample, "not found in window data")
}

# Get founder matrix and sample frequencies
founder_matrix <- window_snps %>% select(all_of(founders))
sample_freqs <- window_snps[[test_sample]]

cat("=== Matrix Analysis ===\n")
cat("Founder matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
cat("Sample frequencies length:", length(sample_freqs), "\n")

# Filter for non-NA values
valid_positions <- !is.na(sample_freqs)
sample_freqs <- sample_freqs[valid_positions]
founder_matrix <- founder_matrix[valid_positions, ]

cat("After filtering NAs:\n")
cat("Valid positions:", sum(valid_positions), "/", length(valid_positions), "\n")
cat("Founder matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")

if (nrow(founder_matrix) > 0) {
  # Convert to matrix for clustering
  founder_matrix <- as.matrix(founder_matrix)
  
  cat("\n=== Testing Sanity Checks ===\n")
  
  # Sanity checks for bad estimation space
  n_snps <- nrow(founder_matrix)
  n_founders <- ncol(founder_matrix)
  
  cat("Number of SNPs:", n_snps, "\n")
  cat("Number of founders:", n_founders, "\n")
  cat("SNPs per founder ratio:", sprintf("%.2f", n_snps / n_founders), "\n\n")
  
  # Check 1: Too few SNPs relative to founders (rule of thumb: need at least 3x)
  cat("Check 1: SNPs per founder ratio\n")
  min_snps_needed <- n_founders * 3
  cat("  Need at least:", min_snps_needed, "SNPs\n")
  cat("  Have:", n_snps, "SNPs\n")
  if (n_snps < min_snps_needed) {
    cat("  ❌ FAILED: Too few SNPs for reliable estimation\n")
    cat("  Recommendation: Increase window size or skip this case\n")
  } else {
    cat("  ✓ PASSED: Sufficient SNPs for estimation\n")
  }
  cat("\n")
  
  # Check 2: Matrix condition number (numerical stability)
  cat("Check 2: Matrix condition number\n")
  if (n_snps >= n_founders) {
    condition_num <- kappa(founder_matrix)
    cat("  Matrix condition number:", format(condition_num, scientific = TRUE), "\n")
    cat("  Threshold: 1e10\n")
    if (condition_num > 1e10) {
      cat("  ❌ FAILED: Matrix is numerically unstable\n")
      cat("  Recommendation: Increase window size or skip this case\n")
    } else {
      cat("  ✓ PASSED: Matrix is numerically stable\n")
    }
  } else {
    cat("  ⚠️  SKIPPED: Not enough SNPs to compute condition number\n")
  }
  cat("\n")
  
  # Check 3: Effective rank (how many founders are actually distinguishable)
  cat("Check 3: Effective matrix rank\n")
  if (n_snps >= n_founders) {
    svd_result <- svd(founder_matrix)
    effective_rank <- sum(svd_result$d > 1e-6)
    cat("  Effective rank:", effective_rank, "\n")
    cat("  Total founders:", n_founders, "\n")
    cat("  Rank ratio:", sprintf("%.2f", effective_rank / n_founders), "\n")
    cat("  Threshold: 0.7\n")
    if (effective_rank < n_founders * 0.7) {
      cat("  ❌ FAILED: Too many founders are indistinguishable\n")
      cat("  Recommendation: Increase window size or skip this case\n")
    } else {
      cat("  ✓ PASSED: Most founders are distinguishable\n")
    }
  } else {
    cat("  ⚠️  SKIPPED: Not enough SNPs to compute effective rank\n")
  }
  cat("\n")
  
  # Overall assessment
  cat("=== Overall Assessment ===\n")
  checks_passed <- 0
  total_checks <- 0
  
  # Count passed checks
  if (n_snps >= min_snps_needed) {
    checks_passed <- checks_passed + 1
  }
  total_checks <- total_checks + 1
  
  if (n_snps >= n_founders) {
    total_checks <- total_checks + 1
    if (kappa(founder_matrix) <= 1e10) {
      checks_passed <- checks_passed + 1
    }
    
    total_checks <- total_checks + 1
    svd_result <- svd(founder_matrix)
    effective_rank <- sum(svd_result$d > 1e-6)
    if (effective_rank >= n_founders * 0.7) {
      checks_passed <- checks_passed + 1
    }
  }
  
  cat("Checks passed:", checks_passed, "/", total_checks, "\n")
  
  if (checks_passed == total_checks) {
    cat("✓ ALL CHECKS PASSED: Safe to proceed with haplotype estimation\n")
  } else {
    cat("❌ SOME CHECKS FAILED: Consider increasing window size or skipping\n")
    cat("   This case would likely produce invalid haplotype frequencies\n")
  }
  
} else {
  cat("✗ No valid data for analysis\n")
}

cat("\n=== Test Complete ===\n")
