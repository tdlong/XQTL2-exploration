#!/usr/bin/env Rscript

# Simple test script to validate haplotype estimation on 100 positions
# This tests that the scripts work correctly before running the full pipeline

cat("=== TESTING HAPLOTYPE ESTIMATION ON 100 POSITIONS ===\n")

# Load libraries
library(tidyverse)
library(limSolve)

# Load parameters
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")

# Load data
cat("Loading REFALT data...\n")
df <- read.table("process/JUICE/RefAlt.chr2R.txt", header = TRUE)

# Transform to frequencies
cat("Converting to frequencies...\n")
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

# Subset dataset
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Get non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("Non-founder samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Test 100 positions in the middle of the chromosome
test_positions <- seq(10000000, 11000000, by = 10000)[1:100]  # 100 positions
window_size_bp <- 20000  # 20kb

cat("Testing", length(test_positions), "positions from 10M to 11M bp\n")
cat("Window size:", window_size_bp, "bp\n\n")

success_count <- 0
total_tested <- 0
results_list <- list()

for (i in seq_along(test_positions)) {
  test_pos <- test_positions[i]
  total_tested <- total_tested + 1
  
  # Show progress every 10 positions
  if (i %% 10 == 0) {
    cat("Progress:", i, "/", length(test_positions), "positions tested,", success_count, "successful\n")
  }
  
  # Define window
  window_start <- test_pos - window_size_bp
  window_end <- test_pos + window_size_bp
  
  # Get SNPs in window
  window_snps <- df3 %>%
    filter(POS >= window_start & POS <= window_end)
  
  # Test first sample
  sample_name <- non_founder_samples[1]
  
  # Get sample data
  sample_data <- window_snps %>%
    filter(name == sample_name) %>%
    select(POS, freq, N)
  
  if (nrow(window_snps) < 10 || nrow(sample_data) < 5) {
    # Return NA for insufficient data
    result_row <- list(
      chr = "chr2R",
      pos = test_pos,
      sample = sample_name,
      window_size = window_size_bp,
      n_snps = nrow(window_snps)
    )
    
    # Add founder frequencies as named columns (all NA)
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- NA
    }
    
    results_list[[length(results_list) + 1]] <- result_row
    next
  }
  
  # Get founder data
  founder_data <- window_snps %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  if (ncol(founder_data) < length(founders) + 1) {
    # Return NA for missing founders
    result_row <- list(
      chr = "chr2R",
      pos = test_pos,
      sample = sample_name,
      window_size = window_size_bp,
      n_snps = nrow(window_snps)
    )
    
    # Add founder frequencies as named columns (all NA)
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- NA
    }
    
    results_list[[length(results_list) + 1]] <- result_row
    next
  }
  
  # Prepare matrices
  founder_matrix <- founder_data %>%
    select(-POS) %>%
    as.matrix()
  
  complete_rows <- complete.cases(founder_matrix)
  
  if (sum(complete_rows) < 5) {
    # Return NA for insufficient complete rows
    result_row <- list(
      chr = "chr2R",
      pos = test_pos,
      sample = sample_name,
      window_size = window_size_bp,
      n_snps = nrow(window_snps)
    )
    
    # Add founder frequencies as named columns (all NA)
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- NA
    }
    
    results_list[[length(results_list) + 1]] <- result_row
    next
  }
  
  founder_matrix <- founder_matrix[complete_rows, ]
  sample_freqs <- sample_data$freq[match(founder_data$POS[complete_rows], sample_data$POS)]
  
  valid_indices <- !is.na(sample_freqs)
  
  if (sum(valid_indices) < 5) {
    # Return NA for insufficient valid frequencies
    result_row <- list(
      chr = "chr2R",
      pos = test_pos,
      sample = sample_name,
      window_size = window_size_bp,
      n_snps = nrow(window_snps)
    )
    
    # Add founder frequencies as named columns (all NA)
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- NA
    }
    
    results_list[[length(results_list) + 1]] <- result_row
    next
  }
  
  founder_matrix <- founder_matrix[valid_indices, ]
  sample_freqs <- sample_freqs[valid_indices]
  
  # Try LSEI
  tryCatch({
    result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, 
                            E = matrix(1, nrow = 1, ncol = ncol(founder_matrix)), 
                            F = 1, G = diag(ncol(founder_matrix)), H = rep(0, ncol(founder_matrix)))
    
    if (result$IsError == 0) {
      success_count <- success_count + 1
      
      # Store result (same format as real script)
      result_row <- list(
        chr = "chr2R",
        pos = test_pos,
        sample = sample_name,
        window_size = window_size_bp,
        n_snps = nrow(window_snps)
      )
      
      # Add founder frequencies as named columns
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- result$X[i]
      }
      
      results_list[[length(results_list) + 1]] <- result_row
    }
    
  }, error = function(e) {
    # Return NA for any error
    result_row <- list(
      chr = "chr2R",
      pos = test_pos,
      sample = sample_name,
      window_size = window_size_bp,
      n_snps = nrow(window_snps)
    )
    
    # Add founder frequencies as named columns (all NA)
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- NA
    }
    
    results_list[[length(results_list) + 1]] <- result_row
  })
}

# Convert results to data frame
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  cat("\n=== TEST RESULTS ===\n")
  cat("Total positions tested:", total_tested, "\n")
  cat("Total results rows:", nrow(results_df), "\n")
  
  # Check if we have the right number of rows (should be total_tested, not just successful)
  expected_rows <- total_tested
  if (nrow(results_df) == expected_rows) {
    cat("✅ Row count is correct!\n")
  } else {
    cat("❌ Row count mismatch! Expected:", expected_rows, "Got:", nrow(results_df), "\n")
  }
  
  # Calculate success rate: positions that return frequency estimates (not NAs)
  successful_estimates <- results_df %>%
    filter(!is.na(B1)) %>%
    nrow()
  
  success_rate <- successful_estimates / total_tested * 100
  
  cat("Successful haplotype estimates (non-NA):", successful_estimates, "\n")
  cat("Success rate:", round(success_rate, 1), "%\n")
  
  # Show first few results
  cat("\nFirst 3 successful estimates:\n")
  if (nrow(results_df) > 0) {
    print(head(results_df, 3))
  }
  
} else {
  cat("\n❌ No successful haplotype estimates!\n")
}

cat("\n=== TEST COMPLETE ===\n")
