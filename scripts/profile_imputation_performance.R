#!/usr/bin/env Rscript

# Profile Imputation Performance
# Test with small subset of SNPs to verify algorithm works and estimate runtime

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "fixed_500kb"  # Test with largest fixed window first
sample_name <- "GJ_3_1"     # Test with just one sample

cat("=== Profiling Imputation Performance ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n")
cat("Test: Random subset of 100 SNPs\n\n")

# Load haplotype results
if (grepl("^fixed_", estimator)) {
  window_size <- as.numeric(gsub("fixed_|kb", "", estimator)) * 1000
  haplotype_file <- file.path(output_dir, paste0("fixed_window_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file) %>%
    filter(window_size == !!window_size)
  cat("✓ Fixed window results loaded for", window_size, "bp window\n")
} else {
  stop("Only fixed estimators supported for now")
}

# Load SNP data
refalt_file <- file.path(output_dir, paste0("df3.", chr, ".RDS"))
df2 <- read_rds(refalt_file)

# Define euchromatin boundaries
euchromatin_start <- 5398184
euchromatin_end <- 24684540

# Filter SNPs to euchromatin
good_snps <- df2 %>%
  filter(name %in% c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")) %>%
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

# Get haplotype positions for this sample
sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name)

haplotype_positions <- sort(unique(sample_haplotypes$pos))
cat("✓ Haplotype positions for sample:", length(haplotype_positions), "\n")
cat("Haplotype range:", min(haplotype_positions), "-", max(haplotype_positions), "bp\n\n")

# Randomly sample 1000 SNPs, ensuring they're within haplotype boundaries
set.seed(123)  # For reproducible sampling
candidate_snps <- valid_snps %>%
  filter(POS > min(haplotype_positions) + 10000,  # Leave 10kb buffer from leftmost
         POS < max(haplotype_positions) - 10000)   # Leave 10kb buffer from rightmost

if (nrow(candidate_snps) < 1000) {
  cat("⚠️  Only", nrow(candidate_snps), "SNPs available for testing\n")
  test_snps <- candidate_snps$POS
} else {
  test_snps <- sample(candidate_snps$POS, 1000)
}

# Sort the test SNPs (as streaming algorithm expects)
test_snps <- sort(test_snps)

cat("=== Test SNP Subset ===\n")
cat("Randomly sampled SNPs:", length(test_snps), "\n")
cat("SNP range:", min(test_snps), "-", max(test_snps), "bp\n")
cat("First 10 SNPs:", paste(head(test_snps, 10), collapse = ", "), "\n")
cat("Last 10 SNPs:", paste(tail(test_snps, 10), collapse = ", "), "\n\n")

# Test the streaming algorithm with this subset
cat("=== Testing Streaming Algorithm ===\n")
start_time <- Sys.time()

# Simplified streaming algorithm for profiling
interpolate_subset <- function(haplotype_results, snp_positions, founders) {
  cat("  Starting interpolation for", length(snp_positions), "SNP positions\n")
  
  # Filter to euchromatin only
  haplotype_results <- haplotype_results %>%
    filter(pos >= euchromatin_start & pos <= euchromatin_end)
  
  cat("  Haplotype results after filter:", nrow(haplotype_results), "rows\n")
  
  # Convert to wide format once
  haplotype_freqs <- haplotype_results %>%
    pivot_wider(names_from = founder, values_from = freq, values_fill = NA)
  
  # Get unique haplotype positions (sorted)
  haplotype_positions <- sort(unique(haplotype_freqs$pos))
  cat("  Unique haplotype positions:", length(haplotype_positions), "\n")
  
  # Sort SNP positions
  snp_positions <- sort(snp_positions)
  
  # Initialize results
  interpolated_results <- list()
  
  # Initialize interval tracking
  current_left_idx <- 1
  current_right_idx <- 2
  
  # Process each SNP sequentially
  for (i in seq_along(snp_positions)) {
    snp_pos <- snp_positions[i]
    
    if (i %% 10 == 0) {
      cat("    Processing SNP", i, "/", length(snp_positions), "at position", snp_pos, "\n")
    }
    
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
    
    interpolated_results[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
  }
  
  cat("  Total interpolated SNPs:", length(interpolated_results), "\n")
  return(interpolated_results)
}

# Run the test
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
test_results <- interpolate_subset(sample_haplotypes, test_snps, founders)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "secs")

cat("\n=== Performance Results ===\n")
cat("Runtime for 1000 SNPs:", round(runtime, 2), "seconds\n")
cat("Rate:", round(1000/runtime, 2), "SNPs per second\n")
cat("Estimated time for full dataset:", round(259659 * runtime / 1000 / 3600, 1), "hours\n")
cat("Estimated time for 6 samples:", round(259659 * 6 * runtime / 1000 / 3600, 1), "hours\n\n")

cat("=== Test Results ===\n")
cat("Successful interpolations:", length(test_results), "/", length(test_snps), "\n")

if (length(test_results) > 0) {
  cat("First few interpolated SNPs:\n")
  first_snp <- names(test_results)[1]
  first_result <- test_results[[1]]
  cat("  SNP", first_snp, ":\n")
  for (founder in names(first_result)) {
    cat("    ", founder, ":", first_result[[founder]][1], "\n")
  }
  
  # Validate imputations against observed frequencies
  cat("\n=== Validation: Imputed vs Observed ===\n")
  
  # Get observed frequencies for test SNPs
  observed_data <- df2 %>%
    filter(name == sample_name, POS %in% test_snps) %>%
    select(POS, freq) %>%
    rename(observed = freq)
  
  # Create comparison table
  validation_data <- data.frame(
    pos = as.numeric(names(test_results)),
    interpolated = sapply(test_results, function(x) {
      # Calculate imputed frequency as weighted sum of founder frequencies
      founder_freqs <- as.numeric(x[1, ])
      # For now, just use the first founder as a simple metric
      # In practice, you'd use the actual founder states for each SNP
      founder_freqs[1]
    })
  ) %>%
    left_join(observed_data, by = c("pos" = "POS"))
  
  # Calculate validation statistics
  valid_comparisons <- validation_data %>%
    filter(!is.na(observed) & !is.na(interpolated))
  
  if (nrow(valid_comparisons) > 0) {
    cat("SNPs with both imputed and observed frequencies:", nrow(valid_comparisons), "\n")
    
    # Calculate correlation and error metrics
    correlation <- cor(valid_comparisons$interpolated, valid_comparisons$observed, use = "complete.obs")
    mae <- mean(abs(valid_comparisons$interpolated - valid_comparisons$observed), na.rm = TRUE)
    rmse <- sqrt(mean((valid_comparisons$interpolated - valid_comparisons$observed)^2, na.rm = TRUE))
    
    cat("Correlation (imputed vs observed):", round(correlation, 3), "\n")
    cat("Mean Absolute Error:", round(mae, 3), "\n")
    cat("Root Mean Square Error:", round(rmse, 3), "\n")
    
    # Show some examples
    cat("\nSample comparisons (first 10):\n")
    print(head(valid_comparisons, 10))
    
    # Check for extreme errors
    extreme_errors <- valid_comparisons %>%
      mutate(error = abs(interpolated - observed)) %>%
      filter(error > 0.5)
    
    if (nrow(extreme_errors) > 0) {
      cat("\n⚠️  SNPs with large errors (>0.5):\n")
      print(head(extreme_errors, 5))
    } else {
      cat("\n✓ No extreme errors found\n")
    }
  } else {
    cat("⚠️  No SNPs with both imputed and observed frequencies for comparison\n")
  }
}

cat("\n=== Recommendations ===\n")
if (runtime < 60) {
  cat("✓ Algorithm works and is reasonably fast\n")
  cat("  Full dataset will take ~", round(259659 * 6 * runtime / 1000 / 3600, 1), "hours\n")
} else if (runtime < 300) {
  cat("⚠️  Algorithm works but is slow\n")
  cat("  Consider optimization or parallelization\n")
} else {
  cat("❌ Algorithm is too slow for practical use\n")
  cat("  Need major optimization or different approach\n")
}
