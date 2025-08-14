#!/usr/bin/env Rscript

# Euchromatic SNP Imputation - Single Estimator
# This script processes one haplotype estimator at a time for parallel processing
# 
# CRITICAL FIX APPLIED: Founder order mismatch in SNP imputation calculation
# - Founder states from pivot_wider were in wrong order (AB8, B1, B2, ...)
# - Haplotype frequencies were in correct order (B1, B2, B3, ..., AB8)
# - This caused systematic calculation errors and -0.016 correlation
# - Fixed by reordering founder states before calculation
#
# Usage: Rscript euchromatic_SNP_imputation_single.R <chr> <param_file> <output_dir> <estimator>

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript euchromatic_SNP_imputation_single.R <chr> <param_file> <output_dir> <estimator>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]

cat("=== Euchromatic SNP Imputation - Single Estimator ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Estimator:", estimator, "\n\n")

# Load parameter file
cat("Loading parameter file...\n")
source(param_file)
cat("✓ Parameter file loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")

# Define euchromatin boundaries for each chromosome
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

# Get boundaries for this chromosome
if (!chr %in% names(euchromatin_boundaries)) {
  stop("Invalid chromosome: ", chr, ". Valid chromosomes: ", paste(names(euchromatin_boundaries), collapse = ", "))
}

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region for", chr, ":", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load haplotype results for this estimator
cat("Loading haplotype estimation results...\n")
if (grepl("^fixed_", estimator)) {
  # Fixed window estimator
  window_size <- as.numeric(gsub("fixed_", "", gsub("kb", "000", estimator)))
  haplotype_file <- file.path(output_dir, paste0("fixed_window_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file) %>%
    filter(window_size == !!window_size)
  cat("✓ Fixed window results loaded for", window_size, "bp window\n")
} else if (grepl("^adaptive_h", estimator)) {
  # Adaptive window estimator
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file) %>%
    filter(h_cutoff == !!h_cutoff)
  cat("✓ Adaptive window results loaded for h_cutoff =", h_cutoff, "\n")
} else {
  stop("Invalid estimator format. Must be 'fixed_<size>kb' or 'adaptive_h<number>'")
}

cat("Haplotype estimates:", nrow(haplotype_results), "\n\n")

# Load observed SNP data from REFALT files
cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data file not found: ", refalt_file)
}

# Load and process REFALT data
df <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(df), "rows\n")

# Transform REF/ALT counts to frequencies (same as haplotype scripts)
cat("Converting counts to frequencies...\n")
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

cat("✓ Processed REFALT data:", nrow(df2), "rows\n")

# Debug: inspect data structure
cat("Data structure:\n")
cat("Columns:", paste(names(df2), collapse = ", "), "\n")
cat("Column types:", paste(sapply(df2, class), collapse = ", "), "\n")
cat("First few rows:\n")
print(head(df2, 3))
cat("\n")

# Data is already processed - no need to convert counts to frequencies
cat("Data already in correct format - skipping transformation\n")

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

# Get valid SNPs for evaluation (euchromatin only)
valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(df2, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n\n")

# Function to interpolate haplotype frequencies to SNP positions using streaming algorithm
interpolate_haplotype_frequencies <- function(haplotype_results, snp_positions, founders) {
  cat("    DEBUG: Starting interpolation for", length(snp_positions), "SNP positions\n")
  
  # Filter to euchromatin only
  haplotype_results <- haplotype_results %>%
    filter(pos >= euchromatin_start & pos <= euchromatin_end)
  
  cat("    DEBUG: Haplotype results after euchromatin filter:", nrow(haplotype_results), "rows\n")
  
  # Check if we have any valid results
  if (nrow(haplotype_results) == 0) {
    cat("    Warning: No haplotype results in euchromatin region\n")
    return(list())
  }
  
  # Check if all frequencies are NA
  all_na_check <- haplotype_results %>%
    filter(!is.na(freq)) %>%
    nrow()
  
  cat("    DEBUG: Non-NA haplotype frequencies:", all_na_check, "rows\n")
  
  if (all_na_check == 0) {
    cat("    Warning: All haplotype frequencies are NA for this sample\n")
    return(list())
  }
  
  # Convert to wide format once
  haplotype_freqs <- haplotype_results %>%
    pivot_wider(names_from = founder, values_from = freq, values_fill = NA)
  
  cat("    DEBUG: Wide format haplotype data:", nrow(haplotype_freqs), "rows\n")
  cat("    DEBUG: Wide format columns:", paste(names(haplotype_freqs), collapse = ", "), "\n")
  
  # Get unique haplotype positions (sorted)
  haplotype_positions <- sort(unique(haplotype_freqs$pos))
  
  cat("    DEBUG: Unique haplotype positions:", length(haplotype_positions), "\n")
  cat("    DEBUG: First few haplotype positions:", paste(head(haplotype_positions, 5), collapse = ", "), "\n")
  cat("    DEBUG: Last few haplotype positions:", paste(tail(haplotype_positions, 5), collapse = ", "), "\n")
  
  # Sort SNP positions
  snp_positions <- sort(snp_positions)
  
  cat("    DEBUG: First few SNP positions:", paste(head(snp_positions, 5), collapse = ", "), "\n")
  cat("    DEBUG: Last few SNP positions:", paste(tail(snp_positions, 5), collapse = ", "), "\n")
  
  # Initialize results
  interpolated_results <- list()
  
  # Initialize interval tracking
  current_left_idx <- 1
  current_right_idx <- 2
  
  # Process each SNP sequentially
  for (snp_pos in snp_positions) {
    # Find the haplotype interval containing this SNP
    while (current_right_idx <= length(haplotype_positions) && 
           haplotype_positions[current_right_idx] < snp_pos) {
      current_left_idx <- current_right_idx
      current_right_idx <- current_right_idx + 1
    }
    
    # Check if we have a valid interval
    if (current_left_idx > length(haplotype_positions) || 
        current_right_idx > length(haplotype_positions)) {
      next  # SNP is beyond haplotype range
    }
    
    left_pos <- haplotype_positions[current_left_idx]
    right_pos <- haplotype_positions[current_right_idx]
    
    # Get founder columns that exist
    existing_founders <- intersect(founders, names(haplotype_freqs))
    
    # Extract frequencies for left and right positions
    left_freqs <- haplotype_freqs %>% 
      filter(pos == left_pos) %>% 
      select(all_of(existing_founders))
    
    right_freqs <- haplotype_freqs %>% 
      filter(pos == right_pos) %>% 
      select(all_of(existing_founders))
    
    # Convert to numeric vectors
    left_freqs_numeric <- as.numeric(left_freqs[1, ])
    right_freqs_numeric <- as.numeric(right_freqs[1, ])
    
    # If both sides are all NA, skip
    if (all(is.na(left_freqs_numeric)) && all(is.na(right_freqs_numeric))) { 
      next 
    }
    
    # Interpolation weight
    alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
    
    # Vectorized interpolation
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
  }
  
  cat("    DEBUG: Total interpolated SNPs:", length(interpolated_results), "\n")
  return(interpolated_results)
}

# Function to calculate imputed SNP frequency from haplotype frequencies and founder states
calculate_imputed_snp_frequency <- function(haplotype_freqs, founder_states) {
  # founder_states should be a vector of 0s and 1s for each founder
  # haplotype_freqs should be a vector of frequencies for each founder
  # imputed_freq = sum(haplotype_freqs * founder_states)
  imputed_freq <- sum(haplotype_freqs * founder_states)
  return(imputed_freq)
}

# Function to create SNP imputation table for this estimator
create_snp_imputation_table <- function(valid_snps, haplotype_results, founders, estimator) {
  cat("Creating SNP imputation table for", estimator, "...\n")
  
  # Get unique SNP positions and samples
  snp_positions <- valid_snps %>% distinct(CHROM, POS) %>% pull(POS)
  all_samples <- unique(valid_snps$name)
  non_founder_samples <- all_samples[!all_samples %in% founders]
  
  cat("SNP positions:", length(snp_positions), "\n")
  cat("Non-founder samples:", paste(non_founder_samples, collapse = ", "), "\n\n")
  
  # Initialize results list
  results_list <- list()
  
  # Process each sample separately
  for (sample_name in non_founder_samples) {
    cat("Processing sample:", sample_name, "\n")
    
    # Get observed frequencies for this sample
    sample_observed <- valid_snps %>%
      filter(name == sample_name) %>%
      select(CHROM, POS, freq) %>%
      rename(observed = freq)
    
    # Get haplotype results for this sample
    sample_haplotypes <- haplotype_results %>%
      filter(sample == sample_name)
    
    if (nrow(sample_haplotypes) == 0) {
      cat("  Warning: No haplotype results for sample", sample_name, "\n")
      next
    }
    
    # OPTIMIZATION: Use the fast streaming algorithm instead of slow interpolation
    cat("  Running fast streaming algorithm for", length(snp_positions), "SNPs...\n")
    
    # Filter to euchromatin only
    sample_haplotypes <- sample_haplotypes %>%
      filter(pos >= 5398184 & pos <= 24684540)
    
    # Convert to wide format once
    haplotype_freqs <- sample_haplotypes %>%
      pivot_wider(names_from = founder, values_from = freq, values_fill = NA)
    
    # Get unique haplotype positions (sorted)
    haplotype_positions <- sort(unique(haplotype_freqs$pos))
    
    # Sort SNP positions
    snp_positions <- sort(snp_positions)
    
    # Initialize results
    interpolated_haplotypes <- list()
    
    # Initialize interval tracking (this is the key to speed!)
    current_left_idx <- 1
    current_right_idx <- 2
    
    # Process each SNP sequentially with efficient interval tracking
    for (snp_idx in seq_along(snp_positions)) {
      snp_pos <- snp_positions[snp_idx]
      
      # Progress indicator every 1000 SNPs
      if (snp_idx %% 1000 == 0) {
        cat("    Processing SNP", snp_idx, "/", length(snp_positions), "at position", snp_pos, "\n")
      }
      
      # Find the haplotype interval containing this SNP (efficient!)
      while (current_right_idx <= length(haplotype_positions) && 
             haplotype_positions[current_right_idx] < snp_pos) {
        current_left_idx <- current_right_idx
        current_right_idx <- current_right_idx + 1
      }
      
      # Check if we have a valid interval
      if (current_left_idx > length(haplotype_positions) || 
          current_right_idx > length(haplotype_positions)) {
        next  # SNP is beyond haplotype range
      }
      
      left_pos <- haplotype_positions[current_left_idx]
      right_pos <- haplotype_positions[current_right_idx]
      
      # Get founder columns that exist - ensure consistent order (FOUNDER ORDER FIX!)
      founder_order <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
      existing_founders <- intersect(founder_order, names(haplotype_freqs))
      existing_founders <- existing_founders[order(match(existing_founders, founder_order))]
      
      # Extract frequencies for left and right positions
      left_freqs <- haplotype_freqs %>% 
        filter(pos == left_pos) %>% 
        select(all_of(existing_founders))
      
      right_freqs <- haplotype_freqs %>% 
        filter(pos == right_pos) %>% 
        select(all_of(existing_founders))
      
      # Convert to numeric vectors
      left_freqs_numeric <- as.numeric(left_freqs[1, ])
      right_freqs_numeric <- as.numeric(right_freqs[1, ])
      
      # If both sides are all NA, skip
      if (all(is.na(left_freqs_numeric)) && all(is.na(right_freqs_numeric))) { 
        next 
      }
      
      # Interpolation weight
      alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
      
      # Vectorized interpolation
      interpolated_freqs <- rep(NA, length(existing_founders))
      
      for (j in seq_along(existing_founders)) {
        left_val <- left_freqs_numeric[j]
        right_val <- right_freqs_numeric[j]
        
        if (is.na(left_val) && is.na(right_val)) {
          interpolated_freqs[j] <- NA
        } else if (is.na(left_val)) {
          interpolated_freqs[j] <- right_val
        } else if (is.na(right_val)) {
          interpolated_freqs[j] <- left_val
        } else {
          interpolated_freqs[j] <- alpha * left_val + (1 - alpha) * right_val
        }
      }
      
      # Store interpolated haplotype frequencies
      interpolated_haplotypes[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
    }
    
    if (length(interpolated_haplotypes) == 0) {
      cat("  Warning: No interpolated haplotypes for sample", sample_name, "\n")
      next
    }
    
    # OPTIMIZATION: Pre-process founder states once for efficient lookup
    cat("  Pre-processing founder states for efficient lookup...\n")
    founder_states_wide <- valid_snps %>%
      filter(name %in% founders) %>%
      select(CHROM, POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq) %>%
      arrange(POS)
    
    # CRITICAL FIX: Ensure founder columns are in correct order
    founder_order <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
    founder_states_wide <- founder_states_wide %>%
      select(CHROM, POS, all_of(founder_order))
    
    cat("  Processing", length(interpolated_haplotypes), "SNPs with optimized lookup...\n")
    
    # Calculate imputed SNP frequencies with vectorized operations
    total_snps <- length(interpolated_haplotypes)
    for (snp_idx in seq_along(interpolated_haplotypes)) {
      snp_pos <- names(interpolated_haplotypes)[snp_idx]
      snp_pos_num <- as.numeric(snp_pos)
      
      # Progress indicator every 1000 SNPs
      if (snp_idx %% 1000 == 0) {
        cat("    Processing SNP", snp_idx, "/", total_snps, "at position", snp_pos, "\n")
      }
      
      # Get observed frequency
      observed_freq <- sample_observed %>%
        filter(POS == snp_pos_num) %>%
        pull(observed)
      
      if (length(observed_freq) == 0) next
      
      # Get haplotype frequencies for this SNP
      haplotype_freqs <- interpolated_haplotypes[[snp_pos]]
      
      # OPTIMIZATION: Efficient founder state lookup using pre-processed table
      snp_founder_states <- founder_states_wide %>%
        filter(POS == snp_pos_num) %>%
        select(all_of(founder_order)) %>%
        as.numeric()
      
      if (length(snp_founder_states) == 0) next
      
      # Founder states are now already in correct order (B1, B2, B3, ..., AB8)
      # No need for reordering - the pivot_wider + select ensures correct order
      
      # Calculate imputed frequency with correctly ordered founder states
      imputed_freq <- calculate_imputed_snp_frequency(haplotype_freqs, snp_founder_states)
      
      # Store result
      results_list[[length(results_list) + 1]] <- list(
        chr = chr,
        pos = snp_pos_num,
        sample = sample_name,
        observed = observed_freq,
        imputed = imputed_freq,
        estimator = estimator
      )
    }
    
    cat("  Completed sample:", sample_name, "\n")
  }
  
  return(results_list)
}

# Create the SNP imputation table
cat("=== Creating SNP Imputation Table ===\n")
start_time <- Sys.time()

imputation_results <- create_snp_imputation_table(valid_snps, haplotype_results, founders, estimator)

if (length(imputation_results) > 0) {
  # Convert to data frame
  imputation_df <- bind_rows(imputation_results)
  
  # Save results
  output_file <- file.path(output_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
  write_rds(imputation_df, output_file)
  
  cat("\n✓ SNP Imputation completed successfully!\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total imputations:", nrow(imputation_df), "\n")
  
  # Summary statistics
  cat("\nSummary by sample:\n")
  summary_stats <- imputation_df %>%
    group_by(sample) %>%
    summarize(
      n_imputations = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_imputed = mean(imputed, na.rm = TRUE),
      correlation = cor(observed, imputed, use = "complete.obs")
    )
  print(summary_stats)
  
} else {
  cat("\n❌ No SNP imputation results obtained!\n")
}

# Timing
total_time <- difftime(Sys.time(), start_time, units = "secs")
cat("\nTotal processing time:", round(total_time, 2), "seconds\n")
cat("Rate:", round(nrow(imputation_df) / as.numeric(total_time), 1), "imputations/second\n")

cat("\n=== SNP Imputation Complete ===\n")
