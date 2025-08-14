#!/usr/bin/env Rscript

# Simple debug script to test haplotype estimation on just a few positions
# This will help us identify why all estimates are NA without scanning the entire chromosome

cat("=== SIMPLE HAPLOTYPE DEBUG ===\n")

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
# Start at 10M bp and go to 20M bp in 100kb steps
test_positions <- seq(10000000, 20000000, by = 100000)
window_size_bp <- 20000  # 20kb

cat("Testing", length(test_positions), "positions from 10M to 20M bp\n")
cat("Window size:", window_size_bp, "bp\n\n")

success_count <- 0
total_tested <- 0

for (i in seq_along(test_positions)) {
  test_pos <- test_positions[i]
  total_tested <- total_tested + 1
  
  # Show progress every 10 positions
  if (i %% 10 == 0) {
    cat("Progress:", i, "/", length(test_positions), "positions tested,", success_count, "successful\n")
  }
  
  # Only show detailed output for first 3 positions
  show_details <- i <= 3
  
  if (show_details) {
    cat("=== Testing position", test_pos, "===\n")
  }
  
  # Define window
  window_start <- test_pos - window_size_bp
  window_end <- test_pos + window_size_bp
  
  # Get SNPs in window
  window_snps <- df3 %>%
    filter(POS >= window_start & POS <= window_end)
  
  if (show_details) {
    cat("  SNPs in window:", nrow(window_snps), "\n")
  }
  
  # Test first sample
  sample_name <- non_founder_samples[1]
  
  if (show_details) {
    cat("  Testing sample:", sample_name, "\n")
  }
  
  # Get sample data
  sample_data <- window_snps %>%
    filter(name == sample_name) %>%
    select(POS, freq, N)
  
  if (show_details) {
    cat("  Sample SNPs:", nrow(sample_data), "\n")
  }
  
  if (nrow(window_snps) < 10 || nrow(sample_data) < 5) {
    if (show_details) {
      cat("  ❌ Insufficient data - window SNPs:", nrow(window_snps), "sample SNPs:", nrow(sample_data), "\n")
    }
    next
  }
  
  # Get founder data
  founder_data <- window_snps %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  if (show_details) {
    cat("  Founder columns found:", ncol(founder_data) - 1, "of", length(founders), "\n")
  }
  
  if (ncol(founder_data) < length(founders) + 2) {
    if (show_details) {
      cat("  ❌ Missing founders\n")
    }
    next
  }
  
  # Prepare matrices
  founder_matrix <- founder_data %>%
    select(-POS) %>%
    as.matrix()
  
  complete_rows <- complete.cases(founder_matrix)
  
  if (show_details) {
    cat("  Complete founder rows:", sum(complete_rows), "\n")
  }
  
  if (sum(complete_rows) < 5) {
    if (show_details) {
      cat("  ❌ Insufficient complete founder rows\n")
    }
    next
  }
  
  founder_matrix <- founder_matrix[complete_rows, ]
  sample_freqs <- sample_data$freq[match(founder_data$POS[complete_rows], sample_data$POS)]
  
  valid_indices <- !is.na(sample_freqs)
  
  if (show_details) {
    cat("  Valid sample frequencies:", sum(valid_indices), "\n")
  }
  
  if (sum(valid_indices) < 5) {
    if (show_details) {
      cat("  ❌ Insufficient valid sample frequencies\n")
    }
    next
  }
  
  founder_matrix <- founder_matrix[valid_indices, ]
  sample_freqs <- sample_freqs[valid_indices]
  
  if (show_details) {
    cat("  Final matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
  }
  
  # Try LSEI
  tryCatch({
    result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, 
                            E = matrix(1, nrow = 1, ncol = ncol(founder_matrix)), 
                            F = 1, G = diag(ncol(founder_matrix)), H = rep(0, ncol(founder_matrix)))
    
    if (result$IsError == 0) {
      success_count <- success_count + 1
      if (show_details) {
        cat("  ✅ SUCCESS! Founder frequencies:", paste(round(result$X, 3), collapse = ", "), "\n")
      }
    } else {
      if (show_details) {
        cat("  ❌ LSEI failed with error code:", result$IsError, "\n")
      }
    }
    
  }, error = function(e) {
    if (show_details) {
      cat("  ❌ LSEI crashed:", e$message, "\n")
    }
  })
  
  if (show_details) {
    cat("\n")
  }
}

cat("\n=== FINAL SUMMARY ===\n")
cat("Total positions tested:", total_tested, "\n")
cat("Successful haplotype estimates:", success_count, "\n")
cat("Success rate:", round(success_count / total_tested * 100, 1), "%\n")
cat("=== DEBUG COMPLETE ===\n")
