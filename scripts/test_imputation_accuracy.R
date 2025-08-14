#!/usr/bin/env Rscript

# =============================================================================
# Test SNP Imputation Accuracy - Simplified Version
# =============================================================================
# This script tests if the founder order fix actually improves accuracy
# Focuses only on the core question: does correlation improve from -0.016?

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
estimator <- "fixed_500kb"
sample_name <- "GJ_3_1"
output_dir <- "process/JUICE"

cat("=== Testing SNP Imputation Accuracy ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n\n")

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

# Get sample haplotypes
sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name)

cat("Sample haplotypes:", nrow(sample_haplotypes), "rows\n")

# Get valid SNPs for this sample
valid_snps <- df2 %>%
  filter(name == sample_name, 
         POS >= euchromatin_start, 
         POS <= euchromatin_end) %>%
  select(POS, freq, N) %>%
  rename(observed_freq = freq, total_read_depth = N)

cat("Valid SNPs for sample:", nrow(valid_snps), "\n\n")

# Test with 1000 SNPs (like the original profiling script)
set.seed(123)
test_snps <- sample(valid_snps$POS, 1000)
test_snps <- sort(test_snps)

cat("=== Testing with 1000 SNPs ===\n")
cat("SNP range:", min(test_snps), "-", max(test_snps), "bp\n\n")

# COPY THE EXACT FAST ALGORITHM from profile_imputation_performance.R
cat("Running the exact fast algorithm from profiling script...\n")
start_time <- Sys.time()

# Filter to euchromatin only
sample_haplotypes <- sample_haplotypes %>%
  filter(pos >= euchromatin_start & pos <= euchromatin_end)

cat("  Haplotype results after filter:", nrow(sample_haplotypes), "rows\n")

# Convert to wide format once
haplotype_freqs <- sample_haplotypes %>%
  pivot_wider(names_from = founder, values_from = freq, values_fill = NA)

# Get unique haplotype positions (sorted)
haplotype_positions <- sort(unique(haplotype_freqs$pos))
cat("  Unique haplotype positions:", length(haplotype_positions), "\n")

# Sort SNP positions
test_snps <- sort(test_snps)

# Initialize results
interpolated_results <- list()

# Initialize interval tracking (this is the key to speed!)
current_left_idx <- 1
current_right_idx <- 2

# Process each SNP sequentially with efficient interval tracking
for (i in seq_along(test_snps)) {
  snp_pos <- test_snps[i]
  
  if (i %% 100 == 0) {
    cat("    Processing SNP", i, "/", length(test_snps), "at position", snp_pos, "\n")
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
  
  # Store interpolated haplotype frequencies (not SNP frequencies yet)
  interpolated_results[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
}

# Now calculate SNP frequencies using the interpolated haplotypes
cat("  Calculating SNP frequencies from interpolated haplotypes...\n")

# Pre-process founder states once for efficiency
founder_states_wide <- df2 %>%
  filter(name %in% c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq) %>%
  arrange(POS) %>%
  select(POS, B1, B2, B3, B4, B5, B6, B7, AB8)

# Calculate imputed SNP frequencies
snp_results <- list()
for (snp_pos in names(interpolated_results)) {
  snp_pos_num <- as.numeric(snp_pos)
  
  # Get founder states for this SNP
  founder_states <- founder_states_wide %>%
    filter(POS == snp_pos_num) %>%
    select(-POS) %>%
    as.numeric()
  
  if (length(founder_states) == 8 && !any(is.na(founder_states))) {
    # Get interpolated haplotype frequencies
    haplotype_freqs <- interpolated_results[[snp_pos]]
    
    # Calculate imputed frequency: sum(haplotype_freq × founder_state)
    imputed_freq <- sum(as.numeric(haplotype_freqs[1, ]) * founder_states)
    
    snp_results[[snp_pos]] <- list(
      pos = snp_pos_num,
      imputed_freq = imputed_freq
    )
  }
}

# Replace interpolated_results with SNP results
interpolated_results <- snp_results

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "secs")

cat("\n=== Performance Results ===\n")
cat("Runtime for 1000 SNPs:", round(as.numeric(runtime), 2), "seconds\n")
cat("Rate:", round(1000/as.numeric(runtime), 2), "SNPs per second\n")
cat("Successfully interpolated SNPs:", length(interpolated_results), "/", length(test_snps), "\n\n")
cat("Estimated time for full dataset (259,659 SNPs):", round(259659 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n")
cat("Estimated time for 6 samples:", round(259659 * 6 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n\n")

# Create simple comparison table
if (length(interpolated_results) > 0) {
  comparison_data <- data.frame(
    pos = sapply(interpolated_results, function(x) x$pos),
    imputed_freq = sapply(interpolated_results, function(x) x$imputed_freq),
    stringsAsFactors = FALSE
  )
  
  # Add observed frequencies
  comparison_data <- comparison_data %>%
    left_join(valid_snps, by = c("pos" = "POS"))
  
  # Calculate correlation
  valid_comparison <- comparison_data %>%
    filter(!is.na(imputed_freq) & !is.na(observed_freq))
  
  if (nrow(valid_comparison) > 10) {
    correlation <- cor(valid_comparison$imputed_freq, valid_comparison$observed_freq, use = "complete.obs")
    cat("=== Accuracy Test ===\n")
    cat("SNPs with both imputed and observed:", nrow(valid_comparison), "\n")
    cat("Correlation (imputed vs observed):", round(correlation, 3), "\n")
    
    if (correlation > 0.5) {
      cat("✅ SUCCESS! Founder order fix worked - correlation improved significantly\n")
    } else if (correlation > 0.1) {
      cat("⚠️  Partial improvement - correlation better but still low\n")
    } else {
      cat("❌ FAILURE! Founder order fix didn't solve the problem\n")
    }
    
    # Show first few results
    cat("\nFirst 5 results:\n")
    print(head(valid_comparison, 5))
    
    # Save full comparison table
    output_file <- paste0("imputation_test_results_", estimator, "_", sample_name, ".csv")
    write_csv(valid_comparison, output_file)
    cat("\n✅ Full results saved to:", output_file, "\n")
    
    # Show summary statistics
    cat("\n=== Summary Statistics ===\n")
    summary_stats <- valid_comparison %>%
      summarise(
        n_snps = n(),
        mean_observed = mean(observed_freq, na.rm = TRUE),
        mean_imputed = mean(imputed_freq, na.rm = TRUE),
        mean_error = mean(abs(imputed_freq - observed_freq), na.rm = TRUE),
        max_error = max(abs(imputed_freq - observed_freq), na.rm = TRUE),
        correlation = correlation
      )
    print(summary_stats)
    
  } else {
    cat("❌ Not enough data for correlation calculation\n")
  }
} else {
  cat("❌ No successful interpolations\n")
}

cat("\n=== Test Complete ===\n")
