#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/create_summary_file_chunked.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== CREATING SUMMARY FILE (SIMPLE APPROACH) ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Define methods
fixed_sizes <- c(20, 50, 100, 200, 500)
h_cutoffs <- c(4, 6, 8, 10)

all_summaries <- list()

# Process fixed window methods
for (size in fixed_sizes) {
  method <- paste0("fixed_", size, "kb")
  cat("Processing", method, "...\n")
  
  # 1. Get SNP counts from haplotype file
  haplo_file <- file.path(results_dir, paste0("fixed_window_", size, "kb_results_", chr, ".RDS"))
  if (!file.exists(haplo_file)) {
    cat("  ❌ Missing haplotype file:", haplo_file, "\n")
    next
  }
  
  haplo_data <- readRDS(haplo_file) %>%
    select(chr, pos, estimate_OK, B1, sample, n_snps) %>%
    mutate(method = method)
  
  # 2. Get MAE from imputation file (average within ±5kb of each position)
  snp_file <- file.path(results_dir, paste0("snp_imputation_fixed_", size, "kb_", chr, ".RDS"))
  if (!file.exists(snp_file)) {
    cat("  ❌ Missing SNP file:", snp_file, "\n")
    next
  }
  
  snp_data <- readRDS(snp_file)
  
  # Calculate MAE from imputation file using proper tidyverse approach
  # Create 10kb bins centered on haplotype positions, then group and summarize
  
  # Get unique haplotype positions for this sample
  haplo_positions <- haplo_data %>%
    select(pos) %>%
    distinct() %>%
    pull(pos)
  
  # Calculate MAE for each haplotype position using proper tidyverse
  mae_data <- snp_data %>%
    # Create bins by finding the closest haplotype position for each SNP
    mutate(
      # Find the closest haplotype position for each SNP
      closest_pos = haplo_positions[sapply(pos, function(x) which.min(abs(haplo_positions - x)))]
    ) %>%
    # Filter to only include SNPs within ±5kb of a haplotype position
    filter(abs(pos - closest_pos) <= 5000) %>%
    # Group by closest position and sample to calculate MAE
    group_by(closest_pos, sample) %>%
    summarize(
      MAE = mean(abs(observed - imputed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Rename for joining
    rename(pos = closest_pos)
  
  # 3. Join everything together
  summary_data <- haplo_data %>%
    left_join(mae_data, by = c("chr", "pos", "sample")) %>%
    mutate(
      NSNPs = n_snps,  # SNP count from haplotype file
      MAE = MAE        # MAE from imputation file
    ) %>%
    select(-n_snps)   # Remove original column
  
  all_summaries[[length(all_summaries) + 1]] <- summary_data
  cat("  ✓ Completed", method, "\n")
}

# Process adaptive window methods
for (h in h_cutoffs) {
  method <- paste0("adaptive_h", h)
  cat("Processing", method, "...\n")
  
  # 1. Get SNP counts from haplotype file
  haplo_file <- file.path(results_dir, paste0("adaptive_window_h", h, "_results_", chr, ".RDS"))
  if (!file.exists(haplo_file)) {
    cat("  ❌ Missing haplotype file:", haplo_file, "\n")
    next
  }
  
  haplo_data <- readRDS(haplo_file) %>%
    select(chr, pos, estimate_OK, B1, sample, n_snps) %>%
    mutate(method = method)
  
  # 2. Get MAE from imputation file (average within ±5kb of each position)
  snp_file <- file.path(results_dir, paste0("snp_imputation_adaptive_h", h, "_", chr, ".RDS"))
  if (!file.exists(snp_file)) {
    cat("  ❌ Missing SNP file:", snp_file, "\n")
    next
  }
  
  snp_data <- readRDS(snp_file)
  
  # Calculate MAE from imputation file using proper tidyverse approach
  # Create 10kb bins centered on haplotype positions, then group and summarize
  
  # Get unique haplotype positions for this sample
  haplo_positions <- haplo_data %>%
    select(pos) %>%
    distinct() %>%
    pull(pos)
  
  # Calculate MAE for each haplotype position using proper tidyverse
  mae_data <- snp_data %>%
    # Create bins by finding the closest haplotype position for each SNP
    mutate(
      # Find the closest haplotype position for each SNP
      closest_pos = haplo_positions[sapply(pos, function(x) which.min(abs(haplo_positions - x)))]
    ) %>%
    # Filter to only include SNPs within ±5kb of a haplotype position
    filter(abs(pos - closest_pos) <= 5000) %>%
    # Group by closest position and sample to calculate MAE
    group_by(closest_pos, sample) %>%
    summarize(
      MAE = mean(abs(observed - imputed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Rename for joining
    rename(pos = closest_pos)
  
  # 3. Join everything together
  summary_data <- haplo_data %>%
    left_join(mae_data, by = c("chr", "pos", "sample")) %>%
    mutate(
      NSNPs = n_snps,  # SNP count from haplotype file
      MAE = MAE        # MAE from imputation file
    ) %>%
    select(-n_snps)   # Remove original column
  
  all_summaries[[length(all_summaries) + 1]] <- summary_data
  cat("  ✓ Completed", method, "\n")
}

# Combine all summaries
final_summary <- bind_rows(all_summaries) %>%
  select(chr, pos, method, B1_freq = B1, estimate_OK, MAE, NSNPs, sample)

# Save summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
saveRDS(final_summary, summary_file)

cat("\n=== FINAL SUMMARY ===\n")
cat("Total rows:", nrow(final_summary), "\n")
cat("Unique positions:", length(unique(final_summary$pos)), "\n")
cat("Unique methods:", length(unique(final_summary$method)), "\n")
cat("Unique samples:", length(unique(final_summary$sample)), "\n")

cat("\n✓ Summary file saved to:", summary_file, "\n")
cat("File size:", round(file.size(summary_file) / 1024^2, 2), "MB\n")
