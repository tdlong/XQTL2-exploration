#!/usr/bin/env Rscript

# =============================================================================
# Test Production SNP Imputation Performance
# =============================================================================
# This script tests if the production code is now fast and accurate

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters (same as production)
chr <- "chr2R"
estimator <- "fixed_500kb"
output_dir <- "process/JUICE"

cat("=== Testing Production SNP Imputation Performance ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load haplotype results
cat("Loading haplotype results...\n")
if (grepl("^fixed_", estimator)) {
  window_size <- as.numeric(gsub("fixed_", "", gsub("kb", "000", estimator)))
  haplotype_file <- file.path(output_dir, paste0("fixed_window_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file) %>%
    filter(window_size == !!window_size)
  cat("✓ Fixed window results loaded for", window_size, "bp window\n")
} else {
  stop("Only testing fixed window for now")
}

cat("Haplotype estimates:", nrow(haplotype_results), "\n\n")

# Load observed SNP data
cat("Loading observed SNP data...\n")
refalt_file <- file.path(output_dir, paste0("df3.", chr, ".RDS"))
df2 <- read_rds(refalt_file)
cat("✓ REFALT data loaded:", nrow(df2), "rows\n")

# Define euchromatin region
euchromatin_start <- 5398184
euchromatin_end <- 24684540

# Filter for high-quality SNPs (same as production)
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

valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end)

cat("✓ Valid euchromatic SNPs:", nrow(valid_snps), "\n")

# Get unique SNP positions and samples
snp_positions <- valid_snps %>% distinct(CHROM, POS) %>% pull(POS)
all_samples <- unique(df2$name)
# Load founders from parameter file
source("helpfiles/JUICE_haplotype_parameters.R")
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("SNP positions:", length(snp_positions), "\n")
cat("Non-founder samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Test with just one sample and 1000 SNPs
sample_name <- "GJ_3_1"
test_snp_positions <- sample(snp_positions, 1000)

cat("=== Testing with sample", sample_name, "and 1000 SNPs ===\n")

# Get observed frequencies for this sample
sample_observed <- df2 %>%
  filter(name == sample_name) %>%
  select(CHROM, POS, freq) %>%
  rename(observed = freq)

# Get haplotype results for this sample
sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name)

if (nrow(sample_haplotypes) == 0) {
  stop("No haplotype results for sample", sample_name)
}

cat("Sample haplotypes:", nrow(sample_haplotypes), "rows\n")

# Test the production algorithm
cat("Running production algorithm...\n")
start_time <- Sys.time()

# OPTIMIZATION: Use the fast streaming algorithm instead of slow interpolation
cat("  Running fast streaming algorithm for", length(test_snp_positions), "SNPs...\n")

# Filter to euchromatin only
sample_haplotypes <- sample_haplotypes %>%
  filter(pos >= euchromatin_start & pos <= euchromatin_end)

# Convert to wide format once
haplotype_freqs <- sample_haplotypes %>%
  pivot_wider(names_from = founder, values_from = freq, values_fill = NA)

# Get unique haplotype positions (sorted)
haplotype_positions <- sort(unique(haplotype_freqs$pos))

# Sort SNP positions
test_snp_positions <- sort(test_snp_positions)

# Initialize results
interpolated_haplotypes <- list()

# Initialize interval tracking (this is the key to speed!)
current_left_idx <- 1
current_right_idx <- 2

# Process each SNP sequentially with efficient interval tracking
for (snp_idx in seq_along(test_snp_positions)) {
  snp_pos <- test_snp_positions[snp_idx]
  
  if (snp_idx %% 100 == 0) {
    cat("    Processing SNP", snp_idx, "/", length(test_snp_positions), "at position", snp_pos, "\n")
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
  stop("No interpolated haplotypes for sample", sample_name)
}

# OPTIMIZATION: Pre-process founder states once for efficient lookup
cat("  Pre-processing founder states for efficient lookup...\n")
founder_states_wide <- df2 %>%
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
results_list <- list()
for (snp_pos in names(interpolated_haplotypes)) {
  snp_pos_num <- as.numeric(snp_pos)
  
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
  imputed_freq <- sum(as.numeric(haplotype_freqs[1, ]) * snp_founder_states)
  
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

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "secs")

cat("\n=== Production Performance Results ===\n")
cat("Runtime for 1000 SNPs:", round(as.numeric(runtime), 2), "seconds\n")
cat("Rate:", round(1000/as.numeric(runtime), 2), "SNPs per second\n")
cat("Successfully interpolated SNPs:", length(results_list), "/", length(test_snp_positions), "\n\n")

# Calculate correlation
if (length(results_list) > 0) {
  # Convert to data frame
  results_df <- bind_rows(results_list)
  
  # Calculate correlation
  correlation <- cor(results_df$observed, results_df$imputed, use = "complete.obs")
  
  cat("=== Production Accuracy Test ===\n")
  cat("SNPs with both imputed and observed:", nrow(results_df), "\n")
  cat("Correlation (imputed vs observed):", round(correlation, 3), "\n")
  
  if (correlation > 0.5) {
    cat("✅ SUCCESS! Production code is both fast AND accurate\n")
  } else {
    cat("❌ FAILURE! Production code has accuracy issues\n")
  }
  
  # Show summary statistics
  cat("\n=== Summary Statistics ===\n")
  summary_stats <- results_df %>%
    summarise(
      n_snps = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_imputed = mean(imputed, na.rm = TRUE),
      mean_error = mean(abs(imputed - observed), na.rm = TRUE),
      max_error = max(abs(imputed - observed), na.rm = TRUE),
      correlation = correlation
    )
  print(summary_stats)
  
  # Estimated time for full dataset
  cat("\n=== Estimated Performance ===\n")
  cat("Estimated time for full dataset (259,659 SNPs):", round(259659 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n")
  cat("Estimated time for 6 samples:", round(259659 * 6 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n")
  
} else {
  cat("❌ No successful interpolations\n")
}

cat("\n=== Production Test Complete ===\n")
