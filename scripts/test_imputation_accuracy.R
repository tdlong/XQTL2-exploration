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

# Test with a small subset first (100 SNPs)
set.seed(123)
test_snps <- sample(valid_snps$POS, 100)
test_snps <- sort(test_snps)

cat("=== Testing with 100 SNPs ===\n")
cat("SNP range:", min(test_snps), "-", max(test_snps), "bp\n\n")

# Simple streaming algorithm for testing
cat("Running simplified streaming algorithm...\n")
start_time <- Sys.time()

# Convert haplotypes to wide format
haplotype_freqs <- sample_haplotypes %>%
  pivot_wider(names_from = founder, values_from = freq, values_fill = NA)

# Get unique haplotype positions
haplotype_positions <- sort(unique(haplotype_freqs$pos))

# Initialize results
interpolated_results <- list()

# Process each test SNP
for (i in seq_along(test_snps)) {
  snp_pos <- test_snps[i]
  
  if (i %% 10 == 0) {
    cat("  Processing SNP", i, "/", length(test_snps), "at position", snp_pos, "\n")
  }
  
  # Find flanking haplotypes
  left_idx <- max(which(haplotype_positions <= snp_pos))
  right_idx <- min(which(haplotype_positions >= snp_pos))
  
  if (is.infinite(left_idx) || is.infinite(right_idx)) {
    next
  }
  
  left_pos <- haplotype_positions[left_idx]
  right_pos <- haplotype_positions[right_idx]
  
  # Get haplotype frequencies at flanking positions
  left_freqs <- haplotype_freqs %>%
    filter(pos == left_pos) %>%
    select(-c(chr, pos, sample, window_size)) %>%
    as.numeric()
  
  right_freqs <- haplotype_freqs %>%
    filter(pos == right_pos) %>%
    select(-c(chr, pos, sample, window_size)) %>%
    as.numeric()
  
  # Interpolate haplotype frequencies
  alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
  interpolated_freqs <- alpha * left_freqs + (1 - alpha) * right_freqs
  
  # Get founder states for this SNP
  founder_states <- df2 %>%
    filter(POS == snp_pos, name %in% c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")) %>%
    select(name, freq) %>%
    pivot_wider(names_from = name, values_from = freq) %>%
    select(B1, B2, B3, B4, B5, B6, B7, AB8) %>%
    as.numeric()
  
  if (length(founder_states) == 8 && !any(is.na(founder_states))) {
    # Calculate imputed frequency: sum(haplotype_freq × founder_state)
    imputed_freq <- sum(interpolated_freqs * founder_states)
    
    interpolated_results[[as.character(snp_pos)]] <- list(
      pos = snp_pos,
      imputed_freq = imputed_freq
    )
  }
}

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "secs")

cat("\n=== Results ===\n")
cat("Runtime:", round(as.numeric(runtime), 2), "seconds\n")
cat("Successfully interpolated SNPs:", length(interpolated_results), "\n\n")

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
