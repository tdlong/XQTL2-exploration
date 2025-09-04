#!/usr/bin/env Rscript

# Run SNP Imputation with List Format Haplotype Results
# Works with the new hap_list_results directory structure

library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/production/run_snp_imputation_list_format.R <chr> <param_file> <output_dir> <method>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
method <- args[4]

# Validate method
if (!method %in% c("adaptive_h4", "smooth_h4")) {
  stop("Method must be 'adaptive_h4' or 'smooth_h4'")
}

cat("=== SNP IMPUTATION WITH LIST FORMAT ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Method:", method, "\n\n")

# Load parameters
source(param_file)

# Load haplotype results in list format
hap_list_results_dir <- file.path(output_dir, "hap_list_results")
haplotype_file <- file.path(hap_list_results_dir, paste0(method, "_list_format_", chr, ".RDS"))

if (!file.exists(haplotype_file)) {
  stop("Haplotype results file not found: ", haplotype_file)
}

haplotype_data <- readRDS(haplotype_file)
cat("✓ Haplotype data loaded:", nrow(haplotype_data), "positions\n")

# Load observed SNP data
observed_file <- file.path(output_dir, "observed_euchromatic", paste0("observed_euchromatic_", chr, ".RDS"))
if (!file.exists(observed_file)) {
  stop("Observed data file not found: ", observed_file)
}

observed_euchromatic <- readRDS(observed_file)
cat("✓ Observed SNP data loaded:", nrow(observed_euchromatic), "rows\n")

# Create results directory
snp_list_results_dir <- file.path(output_dir, "snp_list_results")
dir.create(snp_list_results_dir, showWarnings = FALSE, recursive = TRUE)

# Function to impute SNPs using list format haplotype data
impute_snps_list_format <- function(snp_pos, sample_name, haplotype_data, observed_data, founders) {
  
  # Find the closest haplotype position
  haplotype_positions <- haplotype_data$pos
  closest_idx <- which.min(abs(haplotype_positions - snp_pos))
  closest_pos <- haplotype_positions[closest_idx]
  
  # Check if we're within ±5kb
  if (abs(closest_pos - snp_pos) > 5000) {
    return(NA)  # Too far from any haplotype estimate
  }
  
  # Get haplotype data for this position
  hap_row <- haplotype_data[closest_idx, ]
  
  # Find the sample index
  sample_idx <- which(hap_row$sample[[1]] == sample_name)
  if (length(sample_idx) == 0) {
    return(NA)  # Sample not found
  }
  
  # Get haplotype frequencies for this sample
  hap_freqs <- hap_row$Haps[[1]][[sample_idx]]
  
  # Check if haplotype estimate is valid
  if (any(is.na(hap_freqs))) {
    return(NA)  # Invalid haplotype estimate
  }
  
  # Get observed founder frequencies at this SNP
  snp_data <- observed_data %>%
    filter(POS == snp_pos & name %in% founders)
  
  if (nrow(snp_data) != length(founders)) {
    return(NA)  # Missing founder data
  }
  
  # Create founder frequency vector
  founder_freqs <- snp_data$freq
  names(founder_freqs) <- snp_data$name
  
  # Ensure founders are in the same order
  founder_freqs <- founder_freqs[founders]
  
  # Impute using weighted average of founder frequencies
  imputed_freq <- sum(hap_freqs * founder_freqs)
  
  return(imputed_freq)
}

# Process all SNPs for all samples
cat("Processing SNP imputation...\n")

# Get all unique SNP positions
snp_positions <- unique(observed_euchromatic$POS)
samples <- names_in_bam

# Create imputation results
imputation_results <- expand_grid(
  pos = snp_positions,
  sample = samples
) %>%
  mutate(
    observed = map2_dbl(pos, sample, function(p, s) {
      obs_data <- observed_euchromatic %>%
        filter(POS == p & name == s)
      if (nrow(obs_data) == 1) {
        return(obs_data$freq)
      } else {
        return(NA)
      }
    }),
    imputed = map2_dbl(pos, sample, function(p, s) {
      impute_snps_list_format(p, s, haplotype_data, observed_euchromatic, founders)
    }),
    estimator = method
  ) %>%
  filter(!is.na(observed))  # Only keep SNPs where we have observed data

cat("✓ SNP imputation complete:", nrow(imputation_results), "SNPs processed\n")

# Save results
output_file <- file.path(snp_list_results_dir, paste0("snp_imputation_", method, "_list_format_", chr, ".RDS"))
saveRDS(imputation_results, output_file)

cat("✓ Results saved to:", output_file, "\n")
cat("\n=== SNP IMPUTATION COMPLETE ===\n")
