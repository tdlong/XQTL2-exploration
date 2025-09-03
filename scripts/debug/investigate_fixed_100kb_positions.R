#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/debug/investigate_fixed_100kb_positions.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== INVESTIGATING FIXED 100KB POSITIONS ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n")
cat("Target positions: 8710000, 8720000, 8730000, 8740000\n\n")

# Target positions to investigate
target_positions <- c(8710000, 8720000, 8730000, 8740000)

# 1. Check haplotype file for fixed 100kb
cat("1. Checking fixed 100kb haplotype file...\n")
haplo_file <- file.path(results_dir, "fixed_window_100kb_results_chr2R.RDS")
if (!file.exists(haplo_file)) {
  cat("  ❌ Missing haplotype file:", haplo_file, "\n")
  stop("Cannot proceed without haplotype file")
}

haplo_data <- readRDS(haplo_file)
cat("  ✓ Loaded haplotype file with", nrow(haplo_data), "rows\n")

# Check target positions in haplotype data
target_haplo <- haplo_data %>%
  filter(pos %in% target_positions) %>%
  select(chr, pos, estimate_OK, B1, sample, n_snps)

cat("\n2. Haplotype data for target positions:\n")
if (nrow(target_haplo) > 0) {
  print(target_haplo)
} else {
  cat("  ❌ No haplotype data found for target positions!\n")
}

# 3. Check SNP imputation file for fixed 100kb
cat("\n3. Checking fixed 100kb SNP imputation file...\n")
snp_file <- file.path(results_dir, "snp_imputation_fixed_100kb_chr2R.RDS")
if (!file.exists(snp_file)) {
  cat("  ❌ Missing SNP file:", snp_file, "\n")
  stop("Cannot proceed without SNP file")
}

snp_data <- readRDS(snp_file)
cat("  ✓ Loaded SNP file with", nrow(snp_data), "rows\n")

# Check if there are SNPs around target positions
cat("\n4. Checking for SNPs around target positions (±5kb):\n")
for (pos in target_positions) {
  pos_snps <- snp_data %>%
    filter(pos >= pos - 5000 & pos <= pos + 5000) %>%
    select(chr, pos, sample, observed, imputed)
  
  cat("  Position", pos, ":", nrow(pos_snps), "SNPs found\n")
  if (nrow(pos_snps) > 0) {
    cat("    Sample range:", paste(range(pos_snps$sample), collapse = " to "), "\n")
    cat("    Position range:", paste(range(pos_snps$pos), collapse = " to "), "\n")
  }
}

# 5. Check haplotype positions around target region
cat("\n5. Checking haplotype positions around target region (±50kb):\n")
region_haplo <- haplo_data %>%
  filter(pos >= 8700000 & pos <= 8750000) %>%
  select(chr, pos, estimate_OK, B1, sample, n_snps) %>%
  arrange(pos)

cat("  Haplotype positions in region:", nrow(region_haplo), "rows\n")
if (nrow(region_haplo) > 0) {
  cat("  Position range:", paste(range(region_haplo$pos), collapse = " to "), "\n")
  cat("  estimate_OK summary:\n")
  print(table(region_haplo$estimate_OK, useNA = "ifany"))
  
  cat("\n  First few positions:\n")
  print(head(region_haplo, 10))
  
  cat("\n  Last few positions:\n")
  print(tail(region_haplo, 10))
}

# 6. Check for any SNPs in the entire region
cat("\n6. Checking for SNPs in entire region (±50kb):\n")
region_snps <- snp_data %>%
  filter(pos >= 8700000 & pos <= 8750000) %>%
  select(chr, pos, sample, observed, imputed)

cat("  SNPs in region:", nrow(region_snps), "rows\n")
if (nrow(region_snps) > 0) {
  cat("  Position range:", paste(range(region_snps$pos), collapse = " to "), "\n")
  cat("  Sample range:", paste(range(region_snps$sample), collapse = " to "), "\n")
  
  # Check unique positions
  unique_pos <- unique(region_snps$pos)
  cat("  Unique SNP positions:", length(unique_pos), "\n")
  cat("  First few SNP positions:", paste(head(sort(unique_pos), 10), collapse = ", "), "\n")
}

# 7. Summary and conclusions
cat("\n=== SUMMARY & CONCLUSIONS ===\n")
cat("Target positions:", paste(target_positions, collapse = ", "), "\n")

if (nrow(target_haplo) > 0) {
  cat("✓ Haplotype estimates exist for target positions\n")
  cat("  - estimate_OK values:", paste(unique(target_haplo$estimate_OK), collapse = ", "), "\n")
  cat("  - n_snps values:", paste(unique(target_haplo$n_snps), collapse = ", "), "\n")
} else {
  cat("❌ No haplotype estimates for target positions\n")
}

if (nrow(region_snps) > 0) {
  cat("✓ SNPs exist in the region for imputation\n")
} else {
  cat("❌ No SNPs in region for imputation\n")
}

cat("\n=== NEXT STEPS ===\n")
cat("1. Verify if haplotype estimation actually failed at these positions\n")
cat("2. Check if SNP imputation is working correctly\n")
cat("3. Investigate why we have haplotype estimates when n_snps suggests failure\n")
cat("4. Fix plotting to not show unreliable haplotype estimates\n")

# 8. Create simple table for first sample (Rep01_W_F)
cat("\n=== SIMPLE TABLE FOR Rep01_W_F ===\n")
first_sample <- "Rep01_W_F"
cat("Position | NSNPs | estimate_OK\n")
cat("---------|-------|------------\n")

# Get the 6 positions: 4 target + 2 flanking
positions_to_check <- c(8700000, 8710000, 8720000, 8730000, 8740000, 8750000)

for (pos in positions_to_check) {
  pos_data <- haplo_data %>%
    filter(pos == !!pos, sample == first_sample) %>%
    select(pos, n_snps, estimate_OK)
  
  if (nrow(pos_data) > 0) {
    cat(sprintf("%8.0f | %5d | %10d\n", 
                pos_data$pos, 
                pos_data$n_snps, 
                pos_data$estimate_OK))
  } else {
    cat(sprintf("%8.0f | %5s | %10s\n", pos, "NA", "NA"))
  }
}

# 9. Analyze SNP imputation in the 5 intervals
cat("\n=== SNP IMPUTATION ANALYSIS FOR 5 INTERVALS ===\n")
cat("Interval | Left_OK | Right_OK | SNPs_in_interval | Can_impute?\n")
cat("---------|---------|----------|------------------|------------\n")

for (i in 1:5) {
  left_pos <- positions_to_check[i]
  right_pos <- positions_to_check[i + 1]
  
  # Get estimate_OK for left and right positions
  left_ok <- haplo_data %>%
    filter(pos == left_pos, sample == first_sample) %>%
    pull(estimate_OK)
  
  right_ok <- haplo_data %>%
    filter(pos == right_pos, sample == first_sample) %>%
    pull(estimate_OK)
  
  # Count SNPs in this interval
  interval_snps <- snp_data %>%
    filter(sample == first_sample) %>%
    filter(pos > left_pos & pos < right_pos) %>%
    nrow()
  
  # Determine if imputation is possible
  can_impute <- ifelse(length(left_ok) > 0 && length(right_ok) > 0,
                       ifelse(left_ok == 1 || right_ok == 1, "YES", "NO"),
                       "UNKNOWN")
  
  cat(sprintf("%8.0f-%8.0f | %7d | %9d | %16d | %10s\n",
              left_pos, right_pos, 
              ifelse(length(left_ok) > 0, left_ok, NA),
              ifelse(length(right_ok) > 0, right_ok, NA),
              interval_snps,
              can_impute))
}

# 10. Check specific SNPs in each interval
cat("\n=== DETAILED SNP ANALYSIS BY INTERVAL ===\n")
for (i in 1:5) {
  left_pos <- positions_to_check[i]
  right_pos <- positions_to_check[i + 1]
  
  cat(sprintf("\nInterval %d: %8.0f - %8.0f\n", i, left_pos, right_pos))
  
  # Get SNPs in this interval
  interval_snps <- snp_data %>%
    filter(sample == first_sample) %>%
    filter(pos > left_pos & pos < right_pos) %>%
    select(pos, observed, imputed) %>%
    arrange(pos)
  
  if (nrow(interval_snps) > 0) {
    cat("  SNPs found:", nrow(interval_snps), "\n")
    cat("  Position range:", paste(range(interval_snps$pos), collapse = " to "), "\n")
    cat("  First few SNPs:\n")
    print(head(interval_snps, 5))
    
    # Check for NA imputed values
    na_count <- sum(is.na(interval_snps$imputed))
    cat("  SNPs with imputed = NA:", na_count, "/", nrow(interval_snps), "\n")
  } else {
    cat("  No SNPs found in interval\n")
  }
}
