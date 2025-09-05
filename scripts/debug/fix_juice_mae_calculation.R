#!/usr/bin/env Rscript

# Fix script for JUICE MAE calculation issue
# This addresses the most likely cause: ±5kb filter being too restrictive

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript fix_juice_mae_calculation.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

cat("=== FIXING JUICE MAE CALCULATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

# Define methods to fix
fixed_sizes <- c(20, 50, 100, 200, 500)
h_cutoffs <- c(4, 6, 8, 10)

all_summaries <- list()

# Process fixed window methods with relaxed distance filter
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
  
  # 2. Get MAE from imputation file with relaxed distance filter
  snp_file <- file.path(results_dir, paste0("snp_imputation_fixed_", size, "kb_", chr, ".RDS"))
  if (!file.exists(snp_file)) {
    cat("  ❌ Missing SNP file:", snp_file, "\n")
    next
  }
  
  snp_data <- readRDS(snp_file)
  
  # Get unique haplotype positions for this sample
  haplo_positions <- haplo_data %>%
    select(pos) %>%
    distinct() %>%
    pull(pos)
  
  # Calculate MAE with relaxed distance filter (±10kb instead of ±5kb)
  mae_data <- snp_data %>%
    mutate(
      closest_pos = haplo_positions[sapply(pos, function(x) which.min(abs(haplo_positions - x)))]
    ) %>%
    filter(abs(pos - closest_pos) <= 10000) %>%  # Increased from 5000 to 10000
    group_by(closest_pos, sample) %>%
    summarize(
      MAE = mean(abs(observed - imputed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(chr = haplo_data$chr[1]) %>%
    rename(pos = closest_pos)
  
  # 3. Join everything together
  summary_data <- haplo_data %>%
    left_join(mae_data, by = c("chr", "pos", "sample")) %>%
    mutate(
      NSNPs = ifelse(estimate_OK, n_snps, NA),
      MAE = MAE
    ) %>%
    select(-n_snps)
  
  all_summaries[[length(all_summaries) + 1]] <- summary_data
  cat("  ✓ Completed", method, "\n")
}

# Process adaptive window methods with relaxed distance filter
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
  
  # 2. Get MAE from imputation file with relaxed distance filter
  snp_file <- file.path(results_dir, paste0("snp_imputation_adaptive_h", h, "_", chr, ".RDS"))
  if (!file.exists(snp_file)) {
    cat("  ❌ Missing SNP file:", snp_file, "\n")
    next
  }
  
  snp_data <- readRDS(snp_file)
  
  # Get unique haplotype positions for this sample
  haplo_positions <- haplo_data %>%
    select(pos) %>%
    distinct() %>%
    pull(pos)
  
  # Calculate MAE with relaxed distance filter (±10kb instead of ±5kb)
  mae_data <- snp_data %>%
    mutate(
      closest_pos = haplo_positions[sapply(pos, function(x) which.min(abs(haplo_positions - x)))]
    ) %>%
    filter(abs(pos - closest_pos) <= 10000) %>%  # Increased from 5000 to 10000
    group_by(closest_pos, sample) %>%
    summarize(
      MAE = mean(abs(observed - imputed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(chr = haplo_data$chr[1]) %>%
    rename(pos = closest_pos)
  
  # 3. Join everything together
  summary_data <- haplo_data %>%
    left_join(mae_data, by = c("chr", "pos", "sample")) %>%
    mutate(
      NSNPs = ifelse(estimate_OK, n_snps, NA),
      MAE = MAE
    ) %>%
    select(-n_snps)
  
  all_summaries[[length(all_summaries) + 1]] <- summary_data
  cat("  ✓ Completed", method, "\n")
}

# Combine all summaries
final_summary <- bind_rows(all_summaries) %>%
  select(chr, pos, method, B1_freq = B1, estimate_OK, MAE, NSNPs, sample)

# Save fixed summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, "_fixed.RDS"))
saveRDS(final_summary, summary_file)

cat("\n=== FIXED SUMMARY ===\n")
cat("Total rows:", nrow(final_summary), "\n")
cat("Unique positions:", length(unique(final_summary$pos)), "\n")
cat("Unique methods:", length(unique(final_summary$method)), "\n")
cat("Unique samples:", length(unique(final_summary$sample)), "\n")

# Check MAE coverage
mae_coverage <- final_summary %>%
  group_by(method) %>%
  summarize(
    total_positions = n(),
    positions_with_mae = sum(!is.na(MAE)),
    mae_coverage_pct = round(positions_with_mae / total_positions * 100, 1),
    .groups = "drop"
  )

cat("\nMAE Coverage by Method:\n")
print(mae_coverage)

cat("\n✓ Fixed summary file saved to:", summary_file, "\n")
cat("File size:", round(file.size(summary_file) / 1024^2, 2), "MB\n")

cat("\n=== NEXT STEPS ===\n")
cat("1. Test the fixed summary with the plotting script:\n")
cat("   Rscript scripts/production/plot_summary_region.R", chr, param_file, output_dir, "870\n")
cat("2. If it works, replace the original summary file:\n")
cat("   mv", summary_file, file.path(results_dir, paste0("summary_", chr, ".RDS")), "\n")
