#!/usr/bin/env Rscript

# =============================================================================
# Euchromatic SNP Imputation
# =============================================================================
# Uses haplotype estimates to impute SNP frequencies in euchromatic regions
# Focuses only on euchromatic regions (pos 5,398,184 to 24,684,540)
# Output: chr | pos | sample | observed | fixed_10kb | fixed_20kb | ... | adaptive_h2 | adaptive_h4 | ...

library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript euchromatic_SNP_imputation.R <chromosome> <parameter_file> <output_directory>")
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]

cat("=== Euchromatic SNP Imputation ===\n")
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

# Load observed SNP data
cat("Loading observed SNP data...\n")
refalt_file <- paste0(mydir, "/RefAlt.", mychr, ".txt")
if (!file.exists(refalt_file)) {
  stop("REFALT file not found: ", refalt_file)
}

df <- read.table(refalt_file, header = TRUE)
cat("✓ REFALT data loaded:", nrow(df), "rows\n")

# Transform REF/ALT counts to frequencies
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

# Filter for high-quality SNPs (same criteria as haplotype estimation)
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

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n")

# Function to interpolate haplotype frequencies to SNP positions
interpolate_haplotype_frequencies <- function(haplotype_results, snp_positions, founders) {
  # Filter to euchromatin only
  haplotype_results <- haplotype_results %>%
    filter(pos >= euchromatin_start & pos <= euchromatin_end)
  
  # Check if we have any valid results
  if (nrow(haplotype_results) == 0) {
    cat("  Warning: No haplotype results in euchromatin region\n")
    return(list())
  }
  
  # Check if all frequencies are NA
  all_na_check <- haplotype_results %>%
    filter(!is.na(freq)) %>%
    nrow()
  
  if (all_na_check == 0) {
    cat("  Warning: All haplotype frequencies are NA for this sample\n")
    return(list())
  }
  
  # Get unique positions where we have haplotype estimates
  haplotype_positions <- unique(haplotype_results$pos)
  haplotype_positions <- sort(haplotype_positions)
  
  # Initialize results
  interpolated_results <- list()
  
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
        
        # Check if we have numeric data for interpolation
        if (all(is.na(left_freqs)) || all(is.na(right_freqs))) {
          next
        }
        
        # Convert to numeric and check for valid data
        left_freqs_numeric <- as.numeric(left_freqs)
        right_freqs_numeric <- as.numeric(right_freqs)
        
        if (all(is.na(left_freqs_numeric)) || all(is.na(right_freqs_numeric))) {
          next
        }
        
        # Interpolation weight
        alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
        
        # Linear interpolation: H{snp} = α*H{left} + (1-α)*H{right}
        interpolated_freqs <- alpha * left_freqs_numeric + (1 - alpha) * right_freqs_numeric
        interpolated_results[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
      }
    }
  }
  
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

# Function to create wide format SNP imputation table
create_wide_imputation_table <- function(valid_snps, fixed_results, adaptive_results, founders) {
  cat("Creating wide format SNP imputation table...\n")
  
  # Get unique SNP positions and samples
  snp_positions <- valid_snps %>% distinct(CHROM, POS) %>% pull(POS)
  all_samples <- unique(valid_snps$name)
  non_founder_samples <- all_samples[!all_samples %in% founders]
  
  # Get unique window sizes and h_cutoffs
  fixed_window_sizes <- unique(fixed_results$window_size)
  adaptive_h_cutoffs <- unique(adaptive_results$h_cutoff)
  
  cat("Fixed window sizes:", paste(fixed_window_sizes, collapse = ", "), "\n")
  cat("Adaptive h_cutoffs:", paste(adaptive_h_cutoffs, collapse = ", "), "\n")
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
    
    # Filter haplotype results for this sample
    sample_fixed_results <- fixed_results %>% filter(sample == sample_name)
    sample_adaptive_results <- adaptive_results %>% filter(sample == sample_name)
    
    # Interpolate haplotype frequencies for this sample
    fixed_interpolated <- interpolate_haplotype_frequencies(sample_fixed_results, snp_positions, founders)
    adaptive_interpolated <- interpolate_haplotype_frequencies(sample_adaptive_results, snp_positions, founders)
    
    # Process each SNP position
    for (snp_pos in snp_positions) {
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
    }
  }
  
  return(results_list)
}

# Create SNP imputation table
imputation_results <- create_wide_imputation_table(valid_snps, fixed_results, adaptive_results, founders)

# Convert to data frame
if (length(imputation_results) > 0) {
  imputation_df <- bind_rows(imputation_results)
  
  # Save results
  output_file <- paste0(mydir, "/snp_imputation_euchromatin_", mychr, ".RDS")
  saveRDS(imputation_df, output_file)
  
  cat("\n=== SNP IMPUTATION COMPLETE ===\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total SNP/sample combinations:", nrow(imputation_df), "\n")
  cat("Output format: chr | pos | sample | observed | fixed_10kb | ... | adaptive_h2 | ...\n")
  
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
  
} else {
  cat("\nNo SNP imputation results obtained!\n")
}
