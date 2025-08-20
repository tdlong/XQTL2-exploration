#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/debug_position_563.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== DEBUGGING POSITION 563 ===\n")
cat("Chromosome:", chr, "\n")
cat("Position: 5,630,000 (563 in 10kb units)\n\n")

# Load haplotype results for all methods
methods <- c("fixed_20kb", "fixed_50kb", "fixed_100kb", "fixed_200kb", "fixed_500kb", 
             "adaptive_h4", "adaptive_h6", "adaptive_h8", "adaptive_h10")

haplotype_data <- data.frame()
snp_data <- data.frame()

for (method in methods) {
  # Extract estimator name
  if (grepl("fixed_", method)) {
    size <- gsub("fixed_", "", method)
    estimator <- paste0("fixed_", size)
    file_name <- paste0("fixed_window_", size, "kb_results_", chr, ".RDS")
  } else {
    h <- gsub("adaptive_h", "", method)
    estimator <- paste0("adaptive_h", h)
    file_name <- paste0("adaptive_window_h", h, "_results_", chr, ".RDS")
  }
  
  haplo_file <- file.path(results_dir, file_name)
  snp_file <- file.path(results_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
  
  if (file.exists(haplo_file)) {
    haplo_results <- readRDS(haplo_file) %>%
      filter(pos == 5630000) %>%  # Exact position
      mutate(method = method)
    haplotype_data <- rbind(haplotype_data, haplo_results)
    cat("✓ Loaded haplotype data for", method, "\n")
  } else {
    cat("❌ Missing haplotype file:", haplo_file, "\n")
  }
  
  if (file.exists(snp_file)) {
    snp_results <- readRDS(snp_file) %>%
      filter(pos >= 5625000 & pos <= 5635000) %>%  # ±5kb around position
      mutate(method = method)
    snp_data <- rbind(snp_data, snp_results)
    cat("✓ Loaded SNP data for", method, "\n")
  } else {
    cat("❌ Missing SNP file:", snp_file, "\n")
  }
}

# Analyze haplotype estimates at position 563
cat("\n=== HAPLOTYPE ESTIMATES AT POSITION 563 ===\n")
if (nrow(haplotype_data) > 0) {
  haplo_summary <- haplotype_data %>%
    select(method, sample, B1, B2, B3, B4, B5, B6, B7, B8, estimate_OK) %>%
    arrange(method, sample)
  
  print(haplo_summary)
  
  # Check for identical haplotype estimates across methods
  cat("\n=== CHECKING FOR IDENTICAL HAPLOTYPE ESTIMATES ===\n")
  for (sample_name in unique(haplotype_data$sample)) {
    sample_data <- haplotype_data %>% filter(sample == sample_name)
    if (nrow(sample_data) > 1) {
      # Check if B1 estimates are identical across methods
      b1_values <- sample_data$B1
      if (length(unique(b1_values)) == 1) {
        cat("⚠️  WARNING: B1 estimates identical across all methods for sample", sample_name, "\n")
      } else {
        cat("✓ B1 estimates vary across methods for sample", sample_name, "\n")
        cat("  Range:", min(b1_values), "-", max(b1_values), "\n")
      }
    }
  }
} else {
  cat("❌ No haplotype data found for position 563\n")
}

# Analyze SNP imputation around position 563
cat("\n=== SNP IMPUTATION AROUND POSITION 563 (±5kb) ===\n")
if (nrow(snp_data) > 0) {
  # Calculate RMSE by method
  rmse_by_method <- snp_data %>%
    group_by(method) %>%
    summarise(
      rmse = sqrt(mean((observed - imputed)^2, na.rm = TRUE)),
      n_snps = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_imputed = mean(imputed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(method)
  
  print(rmse_by_method)
  
  # Check for identical RMSE across methods
  if (length(unique(rmse_by_method$rmse)) == 1) {
    cat("\n⚠️  WARNING: RMSE is identical across all methods!\n")
    cat("This suggests either:\n")
    cat("1. Data loading bug (same file loaded multiple times)\n")
    cat("2. SNP imputation algorithm bug (not using haplotype differences)\n")
    cat("3. All methods are equally poor (random baseline)\n")
  } else {
    cat("\n✓ RMSE varies across methods\n")
  }
  
  # Check for identical SNP data across methods
  cat("\n=== CHECKING FOR IDENTICAL SNP DATA ===\n")
  for (sample_name in unique(snp_data$sample)) {
    sample_snp_data <- snp_data %>% filter(sample == sample_name)
    if (nrow(sample_snp_data) > 1) {
      # Check if observed values are identical across methods
      observed_values <- sample_snp_data$observed
      if (length(unique(observed_values)) == 1) {
        cat("⚠️  WARNING: Observed SNP values identical across all methods for sample", sample_name, "\n")
      } else {
        cat("✓ Observed SNP values vary across methods for sample", sample_name, "\n")
      }
      
      # Check if imputed values are identical across methods
      imputed_values <- sample_snp_data$imputed
      if (length(unique(imputed_values)) == 1) {
        cat("⚠️  WARNING: Imputed SNP values identical across all methods for sample", sample_name, "\n")
      } else {
        cat("✓ Imputed SNP values vary across methods for sample", sample_name, "\n")
      }
    }
  }
  
  # Show first few SNPs for each method
  cat("\n=== FIRST 10 SNPS BY METHOD ===\n")
  for (method in unique(snp_data$method)) {
    method_data <- snp_data %>% 
      filter(method == !!method) %>%
      arrange(pos) %>%
      head(10)
    
    cat("\nMethod:", method, "\n")
    print(method_data %>% select(pos, sample, observed, imputed))
  }
  
} else {
  cat("❌ No SNP data found around position 563\n")
}

# Save detailed data for further analysis
debug_file <- file.path(results_dir, paste0("debug_position_563_", chr, ".RDS"))
saveRDS(list(
  haplotype_data = haplotype_data,
  snp_data = snp_data,
  position = 5630000
), debug_file)

cat("\n✓ Debug data saved to:", debug_file, "\n")
