#!/usr/bin/env Rscript

# =============================================================================
# SNP Imputation Diagnostic Script
# =============================================================================
# Tests SNP imputation performance with detailed timing
# Processes subset of data to identify bottlenecks
# Usage: Rscript scripts/snp_imputation_diagnostic.R <chromosome> <parameter_file> <output_directory>

library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript snp_imputation_diagnostic.R <chromosome> <parameter_file> <output_directory>")
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]

cat("=== SNP Imputation Diagnostic ===\n")
cat("Chromosome:", mychr, "\n")
cat("Parameter file:", parfile, "\n")
cat("Output directory:", mydir, "\n\n")

# Source parameter file
if (file.exists(parfile)) {
  source(parfile)
  cat("✓ Parameter file loaded\n")
  cat("Founders:", paste(founders, collapse = ", "), "\n\n")
} else {
  stop("Parameter file not found: ", parfile)
}

# Define euchromatin boundaries
euchromatin_start <- 5398184
euchromatin_end <- 24684540

cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load haplotype estimation results
cat("Loading haplotype estimation results...\n")
start_time <- Sys.time()

# Fixed window results
fixed_file <- paste0(mydir, "/fixed_window_results_", mychr, ".RDS")
if (file.exists(fixed_file)) {
  fixed_results <- readRDS(fixed_file)
  cat("✓ Fixed window results loaded:", nrow(fixed_results), "estimates\n")
} else {
  stop("Fixed window results not found: ", fixed_file)
}

# Adaptive window results
adaptive_file <- paste0(mydir, "/adaptive_window_results_", mychr, ".RDS")
if (file.exists(adaptive_file)) {
  adaptive_results <- readRDS(adaptive_file)
  cat("✓ Adaptive window results loaded:", nrow(adaptive_results), "estimates\n")
} else {
  stop("Adaptive window results not found: ", adaptive_file)
}

load_time <- difftime(Sys.time(), start_time, units = "secs")
cat("Data loading time:", round(load_time, 2), "seconds\n\n")

# Load observed SNP data
cat("Loading observed SNP data...\n")
start_time <- Sys.time()

refalt_file <- paste0(mydir, "/RefAlt.", mychr, ".txt")
if (!file.exists(refalt_file)) {
  stop("REFALT file not found: ", refalt_file)
}

df <- read.table(refalt_file, header = TRUE)
cat("✓ REFALT data loaded:", nrow(df), "rows\n")

# Transform REF/ALT counts to frequencies
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

# Filter for high-quality SNPs
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

# Get valid SNPs for evaluation (euchromatin only)
valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(df2, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n")

data_processing_time <- difftime(Sys.time(), start_time, units = "secs")
cat("Data processing time:", round(data_processing_time, 2), "seconds\n\n")

# Function to interpolate haplotype frequencies to SNP positions (with timing)
interpolate_haplotype_frequencies_timed <- function(haplotype_results, snp_positions, founders, sample_name) {
  cat("  Interpolating for sample:", sample_name, "\n")
  start_time <- Sys.time()
  
  # Filter to euchromatin only
  haplotype_results <- haplotype_results %>%
    filter(pos >= euchromatin_start & pos <= euchromatin_end)
  
  # Check if we have any valid results
  if (nrow(haplotype_results) == 0) {
    cat("    Warning: No haplotype results in euchromatin region\n")
    return(list())
  }
  
  # Check if all frequencies are NA
  all_na_check <- haplotype_results %>%
    filter(!is.na(freq)) %>%
    nrow()
  
  if (all_na_check == 0) {
    cat("    Warning: All haplotype frequencies are NA for this sample\n")
    return(list())
  }
  
  # Get unique positions where we have haplotype estimates
  haplotype_positions <- unique(haplotype_results$pos)
  haplotype_positions <- sort(haplotype_positions)
  
  cat("    Haplotype positions:", length(haplotype_positions), "\n")
  cat("    SNP positions to interpolate:", length(snp_positions), "\n")
  
  # Initialize results
  interpolated_results <- list()
  interpolation_count <- 0
  
  for (snp_pos in snp_positions) {
    # Find nearest haplotype estimates (left and right)
    left_pos <- max(haplotype_positions[haplotype_positions <= snp_pos])
    right_pos <- min(haplotype_positions[haplotype_positions >= snp_pos])
    
    # Check if we have valid positions
    if (is.infinite(left_pos) || is.infinite(right_pos)) {
      next
    }
    
    if (left_pos == right_pos) {
      # SNP is exactly at a haplotype estimate position
      haplotype_freqs <- haplotype_results %>%
        filter(pos == snp_pos) %>%
        select(founder, freq) %>%
        pivot_wider(names_from = founder, values_from = freq)
      
      # Only select founder columns that actually exist
      existing_founders <- intersect(founders, names(haplotype_freqs))
      if (length(existing_founders) > 0) {
        haplotype_freqs <- haplotype_freqs %>% select(all_of(existing_founders))
        interpolated_results[[as.character(snp_pos)]] <- haplotype_freqs
        interpolation_count <- interpolation_count + 1
      }
    } else {
      # Linear interpolation
      left_freqs <- haplotype_results %>%
        filter(pos == left_pos) %>%
        select(founder, freq) %>%
        pivot_wider(names_from = founder, values_from = freq)
      
      right_freqs <- haplotype_results %>%
        filter(pos == right_pos) %>%
        select(founder, freq) %>%
        pivot_wider(names_from = founder, values_from = freq)
      
      # Only select founder columns that actually exist
      existing_founders <- intersect(founders, names(left_freqs))
      if (length(existing_founders) > 0) {
        left_freqs <- left_freqs %>% select(all_of(existing_founders))
        right_freqs <- right_freqs %>% select(all_of(existing_founders))
        
        # Convert to numeric and check for valid data
        left_freqs_numeric <- as.numeric(left_freqs)
        right_freqs_numeric <- as.numeric(right_freqs)
        
        # If either side is all NA, skip this position
        if (all(is.na(left_freqs_numeric)) || all(is.na(right_freqs_numeric))) {
          next
        }
        
        # Interpolation weight
        alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
        
        # Linear interpolation with NA handling
        interpolated_freqs <- rep(NA, length(existing_founders))
        
        for (i in seq_along(existing_founders)) {
          left_val <- left_freqs_numeric[i]
          right_val <- right_freqs_numeric[i]
          
          if (is.na(left_val) && is.na(right_val)) {
            interpolated_freqs[i] <- NA
          } else if (is.na(left_val)) {
            interpolated_freqs[i] <- right_val
          } else if (is.na(right_val)) {
            interpolated_freqs[i] <- left_val
          } else {
            interpolated_freqs[i] <- alpha * left_val + (1 - alpha) * right_val
          }
        }
        
        interpolated_results[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
        interpolation_count <- interpolation_count + 1
      }
    }
  }
  
  interpolation_time <- difftime(Sys.time(), start_time, units = "secs")
  cat("    Interpolation completed:", interpolation_count, "positions in", round(interpolation_time, 2), "seconds\n")
  cat("    Rate:", round(interpolation_count / as.numeric(interpolation_time), 1), "positions/second\n")
  
  return(interpolated_results)
}

# Function to calculate imputed SNP frequency
calculate_imputed_snp_frequency <- function(haplotype_freqs, founder_states) {
  imputed_freq <- sum(haplotype_freqs * founder_states)
  return(imputed_freq)
}

# Test with subset of data
cat("=== Performance Testing ===\n")

# Get subset of samples and positions for testing
all_samples <- unique(valid_snps$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

# Test with first 2 samples and SNP positions that have haplotype estimates
test_samples <- non_founder_samples[1:2]

# Get haplotype positions and find SNPs that are within range
haplotype_positions <- fixed_results %>%
  filter(sample == test_samples[1], pos >= euchromatin_start, pos <= euchromatin_end) %>%
  pull(pos) %>%
  unique() %>%
  sort()

# Find SNPs that are within the haplotype range
test_snp_positions <- valid_snps %>%
  distinct(CHROM, POS) %>%
  filter(POS >= min(haplotype_positions), POS <= max(haplotype_positions)) %>%
  pull(POS) %>%
  head(10)

cat("Testing with", length(test_samples), "samples and", length(test_snp_positions), "SNP positions\n")
cat("Haplotype position range:", min(haplotype_positions), "-", max(haplotype_positions), "bp\n")
cat("SNP position range:", min(test_snp_positions), "-", max(test_snp_positions), "bp\n\n")

# Get unique window sizes and h_cutoffs
fixed_window_sizes <- unique(fixed_results$window_size)
adaptive_h_cutoffs <- unique(adaptive_results$h_cutoff)

cat("Fixed window sizes:", paste(fixed_window_sizes, collapse = ", "), "\n")
cat("Adaptive h_cutoffs:", paste(adaptive_h_cutoffs, collapse = ", "), "\n\n")

# Initialize results list
results_list <- list()

# Process each test sample separately
for (sample_name in test_samples) {
  cat("Processing sample:", sample_name, "\n")
  sample_start_time <- Sys.time()
  
  # Get observed frequencies for this sample
  sample_observed <- valid_snps %>%
    filter(name == sample_name, POS %in% test_snp_positions) %>%
    select(CHROM, POS, freq) %>%
    rename(observed = freq)
  
  # Filter haplotype results for this sample
  sample_fixed_results <- fixed_results %>% filter(sample == sample_name)
  sample_adaptive_results <- adaptive_results %>% filter(sample == sample_name)
  
  # Interpolate haplotype frequencies for this sample
  fixed_interpolated <- interpolate_haplotype_frequencies_timed(sample_fixed_results, test_snp_positions, founders, paste0(sample_name, " (fixed)"))
  adaptive_interpolated <- interpolate_haplotype_frequencies_timed(sample_adaptive_results, test_snp_positions, founders, paste0(sample_name, " (adaptive)"))
  
  # Process each SNP position
  snp_processing_start <- Sys.time()
  snp_count <- 0
  
  for (snp_pos in test_snp_positions) {
    # Start with observed frequency
    observed_freq <- sample_observed %>%
      filter(POS == snp_pos) %>%
      pull(observed)
    
    if (length(observed_freq) == 0) next
    
    # Get founder states for this SNP
    snp_founder_data <- valid_snps %>%
      filter(POS == snp_pos, name %in% founders) %>%
      select(name, freq) %>%
      arrange(match(name, founders))
    
    if (nrow(snp_founder_data) != length(founders)) next
    
    founder_states <- snp_founder_data$freq
    
    # Initialize result row
    result_row <- list(
      chr = mychr,
      pos = snp_pos,
      sample = sample_name,
      observed = observed_freq
    )
    
    # Add fixed window imputed frequencies
    for (window_size in fixed_window_sizes) {
      if (as.character(snp_pos) %in% names(fixed_interpolated)) {
        haplotype_freqs <- fixed_interpolated[[as.character(snp_pos)]]
        if (ncol(haplotype_freqs) == length(founders)) {
          imputed_freq <- calculate_imputed_snp_frequency(as.numeric(haplotype_freqs), founder_states)
          result_row[[paste0("fixed_", window_size/1000, "kb")]] <- imputed_freq
        } else {
          result_row[[paste0("fixed_", window_size/1000, "kb")]] <- NA
        }
      } else {
        result_row[[paste0("fixed_", window_size/1000, "kb")]] <- NA
      }
    }
    
    # Add adaptive window imputed frequencies
    for (h_cutoff in adaptive_h_cutoffs) {
      if (as.character(snp_pos) %in% names(adaptive_interpolated)) {
        haplotype_freqs <- adaptive_interpolated[[as.character(snp_pos)]]
        if (ncol(haplotype_freqs) == length(founders)) {
          imputed_freq <- calculate_imputed_snp_frequency(as.numeric(haplotype_freqs), founder_states)
          result_row[[paste0("adaptive_h", h_cutoff)]] <- imputed_freq
        } else {
          result_row[[paste0("adaptive_h", h_cutoff)]] <- NA
        }
      } else {
        result_row[[paste0("adaptive_h", h_cutoff)]] <- NA
      }
    }
    
    results_list[[length(results_list) + 1]] <- result_row
    snp_count <- snp_count + 1
  }
  
  sample_time <- difftime(Sys.time(), sample_start_time, units = "secs")
  snp_processing_time <- difftime(Sys.time(), snp_processing_start, units = "secs")
  
  cat("  Sample completed:", snp_count, "SNPs in", round(sample_time, 2), "seconds\n")
  cat("  SNP processing time:", round(snp_processing_time, 2), "seconds\n")
  cat("  Rate:", round(snp_count / sample_time, 1), "SNPs/second\n\n")
}

# Convert to data frame
if (length(results_list) > 0) {
  imputation_df <- bind_rows(results_list)
  
  # Save results
  output_file <- paste0(mydir, "/snp_imputation_diagnostic_", mychr, ".RDS")
  saveRDS(imputation_df, output_file)
  
  cat("=== DIAGNOSTIC COMPLETE ===\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total SNP/sample combinations:", nrow(imputation_df), "\n")
  
  # Show summary
  cat("\nSummary by imputation method:\n")
  summary_stats <- imputation_df %>%
    select(-c(chr, pos, sample, observed)) %>%
    summarize(
      across(everything(), ~sum(!is.na(.x)), .names = "n_valid_{.col}"),
      across(everything(), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      across(everything(), ~sd(.x, na.rm = TRUE), .names = "sd_{.col}")
    ) %>%
    pivot_longer(everything(), names_to = c("stat", "method"), names_pattern = "(n_valid|mean|sd)_(.+)") %>%
    pivot_wider(names_from = stat, values_from = value)
  
  print(summary_stats)
  
  # Performance projection
  total_snps <- valid_snps %>% distinct(CHROM, POS) %>% nrow()
  total_samples <- length(non_founder_samples)
  total_combinations <- total_snps * total_samples
  
  cat("\n=== PERFORMANCE PROJECTION ===\n")
  cat("Total SNPs in euchromatin:", total_snps, "\n")
  cat("Total samples:", total_samples, "\n")
  cat("Total combinations:", total_combinations, "\n")
  
  # Estimate time based on test performance
  test_combinations <- nrow(imputation_df)
  test_time <- difftime(Sys.time(), start_time, units = "secs")
  rate_per_combination <- test_time / test_combinations
  estimated_total_time <- rate_per_combination * total_combinations
  
  cat("Test combinations:", test_combinations, "\n")
  cat("Test time:", round(test_time, 2), "seconds\n")
  cat("Rate:", round(rate_per_combination, 4), "seconds per combination\n")
  cat("Estimated total time:", round(estimated_total_time / 3600, 1), "hours\n")
  
} else {
  cat("\nNo SNP imputation results obtained!\n")
}
