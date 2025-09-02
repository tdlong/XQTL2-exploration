#!/usr/bin/env Rscript

# Tidyverse-compliant SNP Imputation Script
# This version follows proper dplyr idioms and is much cleaner

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript euchromatic_SNP_imputation_single_tidy.R <chr> <param_file> <output_dir> <estimator> [start_pos] [max_snps]")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]
test_start_pos <- if (length(args) >= 5) as.numeric(args[5]) else NULL
test_max_snps <- if (length(args) >= 6) as.numeric(args[6]) else NULL

cat("=== Tidyverse SNP Imputation ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Estimator:", estimator, "\n")
if (!is.null(test_start_pos)) cat("Test start position:", test_start_pos, "\n")
if (!is.null(test_max_snps)) cat("Test max SNPs:", test_max_snps, "\n")
cat("\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")

# Define euchromatin boundaries
euchromatin_boundaries <- list(
  chr2R = c(5398184, 24684540)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load haplotype results
cat("Loading haplotype results...\n")
haplotype_file <- if (str_detect(estimator, "^fixed_")) {
  window_size_kb <- str_extract(estimator, "\\d+")
  file.path(output_dir, paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
} else if (str_detect(estimator, "^adaptive_h")) {
  h_cutoff <- str_extract(estimator, "\\d+")
  file.path(output_dir, paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
} else {
  stop("Invalid estimator format")
}

haplotype_results <- read_rds(haplotype_file)
cat("✓ Haplotype results loaded:", nrow(haplotype_results), "rows\n")

# Load REFALT data
cat("Loading REFALT data...\n")
refalt_file <- file.path(dirname(output_dir), paste0("RefAlt.", chr, ".txt"))
refalt_data <- read_tsv(refalt_file, show_col_types = FALSE)
cat("✓ REFALT data loaded:", nrow(refalt_data), "rows\n")

# Process REFALT data using tidyverse
cat("Processing REFALT data...\n")
snp_data <- refalt_data %>%
  # Convert to long format
  pivot_longer(
    cols = -c(CHROM, POS),
    names_to = "sample",
    values_to = "count"
  ) %>%
  # Extract REF/ALT and sample name
  mutate(
    ref_alt = str_sub(sample, 1, 3),
    sample_name = str_sub(sample, 5)
  ) %>%
  select(-sample) %>%
  # Convert to wide format with REF/ALT columns
  pivot_wider(
    names_from = ref_alt,
    values_from = count
  ) %>%
  # Calculate frequencies
  mutate(
    freq = REF / (REF + ALT),
    total_count = REF + ALT
  ) %>%
  select(-c(REF, ALT)) %>%
  # Filter for high-quality SNPs
  group_by(CHROM, POS) %>%
  filter(
    # No missing data in founders
    all(total_count > 0 | !(sample_name %in% founders)),
    # Not fixed in founders
    sum(freq > 0.03 & freq < 0.97, na.rm = TRUE) == 0,
    # Informative
    sum(freq, na.rm = TRUE) > 0.05 | sum(freq, na.rm = TRUE) < 0.95
  ) %>%
  ungroup() %>%
  # Filter to euchromatin
  filter(POS >= euchromatin_start, POS <= euchromatin_end)

cat("✓ Processed SNP data:", n_distinct(snp_data$POS), "SNPs\n")

# Apply test filters if specified
if (!is.null(test_start_pos) || !is.null(test_max_snps)) {
  snp_positions <- snp_data %>%
    distinct(POS) %>%
    pull(POS)
  
  if (!is.null(test_start_pos)) {
    snp_positions <- snp_positions[snp_positions >= test_start_pos]
  }
  
  if (!is.null(test_max_snps) && length(snp_positions) > test_max_snps) {
    snp_positions <- head(sort(snp_positions), test_max_snps)
  }
  
  snp_data <- snp_data %>%
    filter(POS %in% snp_positions)
  
  cat("✓ Test mode: Using", length(snp_positions), "SNPs\n")
}

# Function to interpolate haplotype frequencies using tidyverse
interpolate_haplotype_frequencies <- function(haplotype_data, snp_positions, founder_order) {
  # Filter haplotype data to euchromatin
  haplotype_filtered <- haplotype_data %>%
    filter(pos >= euchromatin_start, pos <= euchromatin_end)
  
  if (nrow(haplotype_filtered) == 0) {
    return(tibble())
  }
  
  # Create interpolation table
  interpolation_data <- tibble(
    snp_pos = snp_positions
  ) %>%
    # Find left and right haplotype positions for each SNP
    mutate(
      left_pos = map_dbl(snp_pos, ~{
        haplotype_filtered %>%
          filter(pos <= .x) %>%
          pull(pos) %>%
          max()
      }),
      right_pos = map_dbl(snp_pos, ~{
        haplotype_filtered %>%
          filter(pos >= .x) %>%
          pull(pos) %>%
          min()
      })
    ) %>%
    # Calculate interpolation weights
    mutate(
      alpha = (right_pos - snp_pos) / (right_pos - left_pos),
      alpha = if_else(is.infinite(alpha) | is.nan(alpha), 0, alpha)
    )
  
  # Get haplotype frequencies for left and right positions
  left_freqs <- haplotype_filtered %>%
    filter(pos %in% interpolation_data$left_pos) %>%
    select(pos, all_of(founder_order)) %>%
    rename_with(~paste0("left_", .x), -pos) %>%
    rename(left_pos = pos)
  
  right_freqs <- haplotype_filtered %>%
    filter(pos %in% interpolation_data$right_pos) %>%
    select(pos, all_of(founder_order)) %>%
    rename_with(~paste0("right_", .x), -pos) %>%
    rename(right_pos = pos)
  
  # Join and interpolate
  interpolated <- interpolation_data %>%
    left_join(left_freqs, by = "left_pos") %>%
    left_join(right_freqs, by = "right_pos") %>%
    # Interpolate each founder frequency
    mutate(across(
      starts_with("left_") & !starts_with("left_pos"),
      ~{
        founder <- str_remove(cur_column(), "left_")
        right_col <- paste0("right_", founder)
        alpha * .x + (1 - alpha) * !!sym(right_col)
      },
      .names = "{.col}"
    )) %>%
    # Select only interpolated frequencies
    select(snp_pos, starts_with("left_") & !starts_with("left_pos")) %>%
    rename_with(~str_remove(.x, "left_"), -snp_pos)
  
  return(interpolated)
}

# Function to calculate imputed SNP frequency
calculate_imputed_frequency <- function(haplotype_freqs, founder_states) {
  sum(haplotype_freqs * founder_states, na.rm = TRUE)
}

# Main imputation function
perform_snp_imputation <- function(snp_data, haplotype_results, founders, estimator) {
  cat("Performing SNP imputation...\n")
  
  # Define founder order
  founder_order <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
  
  # Get non-founder samples
  non_founder_samples <- setdiff(unique(snp_data$sample_name), founders)
  
  # Process each sample
  results <- map_dfr(non_founder_samples, function(sample_name) {
    cat("Processing sample:", sample_name, "\n")
    
    # Get observed frequencies for this sample
    observed_data <- snp_data %>%
      filter(sample_name == !!sample_name) %>%
      select(POS, freq) %>%
      rename(observed = freq)
    
    # Get haplotype results for this sample
    sample_haplotypes <- haplotype_results %>%
      filter(sample == sample_name)
    
    if (nrow(sample_haplotypes) == 0) {
      cat("  Warning: No haplotype results for", sample_name, "\n")
      return(tibble())
    }
    
    # Interpolate haplotype frequencies
    snp_positions <- observed_data$POS
    interpolated_haplotypes <- interpolate_haplotype_frequencies(
      sample_haplotypes, snp_positions, founder_order
    )
    
    if (nrow(interpolated_haplotypes) == 0) {
      cat("  Warning: No interpolated haplotypes for", sample_name, "\n")
      return(tibble())
    }
    
    # Get founder states for each SNP
    founder_states <- snp_data %>%
      filter(sample_name %in% founders) %>%
      select(POS, sample_name, freq) %>%
      pivot_wider(
        names_from = sample_name,
        values_from = freq
      ) %>%
      select(POS, all_of(founder_order))
    
    # Calculate imputed frequencies
    imputation_results <- observed_data %>%
      left_join(interpolated_haplotypes, by = c("POS" = "snp_pos")) %>%
      left_join(founder_states, by = "POS") %>%
      # Calculate imputed frequency for each SNP
      mutate(
        imputed = pmap_dbl(
          list(
            haplotype_freqs = select(cur_data(), all_of(founder_order)),
            founder_states = select(cur_data(), all_of(founder_order))
          ),
          ~calculate_imputed_frequency(as.numeric(..1), as.numeric(..2))
        )
      ) %>%
      select(POS, observed, imputed) %>%
      mutate(
        sample = sample_name,
        estimator = estimator
      )
    
    cat("  Completed:", nrow(imputation_results), "imputations\n")
    return(imputation_results)
  })
  
  return(results)
}

# Run imputation
cat("Starting SNP imputation...\n")
start_time <- Sys.time()

imputation_results <- perform_snp_imputation(snp_data, haplotype_results, founders, estimator)

if (nrow(imputation_results) > 0) {
  # Save results
  output_file <- if (!is.null(test_start_pos) || !is.null(test_max_snps)) {
    file.path(output_dir, paste0("snp_imputation_", estimator, "_", chr, "_TEST.RDS"))
  } else {
    file.path(output_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
  }
  
  write_rds(imputation_results, output_file)
  
  # Summary statistics
  summary_stats <- imputation_results %>%
    group_by(sample) %>%
    summarise(
      n_imputations = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_imputed = mean(imputed, na.rm = TRUE),
      correlation = cor(observed, imputed, use = "complete.obs"),
      rmse = sqrt(mean((observed - imputed)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  cat("\n✓ SNP Imputation completed successfully!\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total imputations:", nrow(imputation_results), "\n\n")
  
  cat("Summary by sample:\n")
  print(summary_stats)
  
} else {
  cat("\n❌ No SNP imputation results obtained!\n")
}

# Timing
total_time <- difftime(Sys.time(), start_time, units = "secs")
cat("\nTotal processing time:", round(total_time, 2), "seconds\n")
cat("Rate:", round(nrow(imputation_results) / as.numeric(total_time), 1), "imputations/second\n")

cat("\n=== Tidyverse SNP Imputation Complete ===\n")
