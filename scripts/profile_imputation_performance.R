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
interpolate_subset <- function(haplotype_results, snp_positions, founders, founder_order) {
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
    
    # Get founder columns that exist - ensure consistent order
existing_founders <- intersect(founders, names(haplotype_freqs))
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
    
    interpolated_results[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
  }
  
  cat("  Total interpolated SNPs:", length(interpolated_results), "\n")
  return(interpolated_results)
}

# Run the test
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
test_results <- interpolate_subset(sample_haplotypes, test_snps, founders, founder_order)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "secs")

# Create comprehensive validation table
cat("\n=== Creating Validation Table ===\n")

# Get founder genotypes and read depths for test SNPs
founder_data <- df2 %>%
  filter(name %in% c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8"),
         POS %in% test_snps) %>%
  select(POS, name, freq, N)

# Define consistent founder order
founder_order <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Check founder order in haplotype data
cat("=== Checking Founder Order Consistency ===\n")
haplotype_cols <- names(sample_haplotypes %>% select(-c(chr, pos, sample, window_size)))
cat("Haplotype data columns:", paste(haplotype_cols, collapse = ", "), "\n")

# Check founder order in SNP data
founder_names_in_data <- unique(founder_data$name)
cat("Founder names in SNP data:", paste(founder_names_in_data, collapse = ", "), "\n")

# Verify order consistency
if (identical(haplotype_cols, founder_order)) {
  cat("✓ Founder order is consistent between haplotype and SNP data\n")
} else {
  cat("⚠️  FOUNDER ORDER MISMATCH!\n")
  cat("  Expected:", paste(founder_order, collapse = ", "), "\n")
  cat("  Haplotype:", paste(haplotype_cols, collapse = ", "), "\n")
  cat("  SNP data:", paste(founder_names_in_data, collapse = ", "), "\n")
}
cat("\n")

# Get observed frequencies and read depths for test SNPs
observed_data <- df2 %>%
  filter(name == sample_name, POS %in% test_snps) %>%
  select(POS, freq, N) %>%
  rename(observed_freq = freq, total_read_depth = N)

# Create validation table with founder states
validation_table <- data.frame(
  chr = chr,
  pos = as.numeric(names(test_results)),
  observed_freq = NA,
  imputed_freq = NA,
  total_read_depth_at_SNP = NA,
  founder_B1_state = NA,
  founder_B2_state = NA,
  founder_B3_state = NA,
  founder_B4_state = NA,
  founder_B5_state = NA,
  founder_B6_state = NA,
  founder_B7_state = NA,
  founder_AB8_state = NA,
  stringsAsFactors = FALSE
)

# Fill in the table
for (i in 1:nrow(validation_table)) {
  snp_pos <- validation_table$pos[i]
  
  # Get observed frequency and read depth
  obs_data <- observed_data %>% filter(POS == snp_pos)
  if (nrow(obs_data) > 0) {
    validation_table$observed_freq[i] <- obs_data$freq
    validation_table$total_read_depth_at_SNP[i] <- obs_data$total_read_depth
  }
  
  # Get founder states for this SNP
  founder_states <- founder_data %>%
    filter(POS == snp_pos) %>%
    select(name, freq)
  
  if (nrow(founder_states) > 0) {
    for (founder in founders) {
      founder_state <- founder_states %>% filter(name == founder) %>% pull(freq)
      if (length(founder_state) > 0) {
        col_name <- paste0("founder_", founder, "_state")
        validation_table[i, col_name] <- founder_state
      }
    }
  }
  
  # Calculate imputed frequency
  if (as.character(snp_pos) %in% names(test_results)) {
    # Get founder genotypes for this SNP
    founder_states <- founder_data %>%
      filter(POS == snp_pos) %>%
      select(name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    # Get interpolated haplotype frequencies
    haplotype_freqs <- test_results[[as.character(snp_pos)]]
    
    if (nrow(haplotype_freqs) > 0 && nrow(founder_states) > 0) {
      # Calculate imputed frequency: sum(haplotype_freq × founder_state)
      imputed_freq <- 0
      for (founder in founders) {
        if (founder %in% names(haplotype_freqs) && founder %in% names(founder_states)) {
          haplo_freq <- haplotype_freqs[[founder]][1]
          founder_state <- founder_states[[founder]][1]
          
          if (!is.na(haplo_freq) && !is.na(founder_state)) {
            imputed_freq <- imputed_freq + (haplo_freq * founder_state)
          }
        }
      }
      validation_table$imputed_freq[i] <- imputed_freq
    }
  }
}

# Calculate errors
validation_table <- validation_table %>%
  mutate(
    error = abs(imputed_freq - observed_freq),
    error_squared = (imputed_freq - observed_freq)^2
  )

# Show summary statistics
cat("Validation table created with", nrow(validation_table), "SNPs\n")
cat("SNPs with both observed and imputed:", sum(!is.na(validation_table$observed_freq) & !is.na(validation_table$imputed_freq)), "\n\n")

# Show error distribution
cat("=== Error Distribution ===\n")
error_summary <- validation_table %>%
  filter(!is.na(error)) %>%
  summarise(
    mean_error = mean(error, na.rm = TRUE),
    median_error = median(error, na.rm = TRUE),
    sd_error = sd(error, na.rm = TRUE),
    q25 = quantile(error, 0.25, na.rm = TRUE),
    q75 = quantile(error, 0.75, na.rm = TRUE),
    max_error = max(error, na.rm = TRUE)
  )
print(error_summary)

# Show correlation by read depth
cat("\n=== Correlation by Read Depth ===\n")
high_coverage <- validation_table %>%
  filter(total_read_depth_at_SNP >= 100, !is.na(imputed_freq), !is.na(observed_freq))

low_coverage <- validation_table %>%
  filter(total_read_depth_at_SNP < 100, !is.na(imputed_freq), !is.na(observed_freq))

if (nrow(high_coverage) > 10) {
  high_corr <- cor(high_coverage$imputed_freq, high_coverage$observed_freq, use = "complete.obs")
  cat("High coverage (≥100 reads) correlation:", round(high_corr, 3), "n =", nrow(high_coverage), "\n")
}

if (nrow(low_coverage) > 10) {
  low_corr <- cor(low_coverage$imputed_freq, low_coverage$observed_freq, use = "complete.obs")
  cat("Low coverage (<100 reads) correlation:", round(low_corr, 3), "n =", nrow(low_coverage), "\n")
}

# Save validation table
output_file <- file.path(output_dir, paste0("validation_table_", estimator, "_", chr, ".csv"))
write_csv(validation_table, output_file)
cat("\nValidation table saved to:", output_file, "\n")

cat("\n=== Performance Results ===\n")
cat("Runtime for 1000 SNPs:", round(runtime, 2), "seconds\n")
cat("Rate:", round(1000/as.numeric(runtime), 2), "SNPs per second\n")
cat("Estimated time for full dataset:", round(259659 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n")
cat("Estimated time for 6 samples:", round(259659 * 6 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n\n")

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
if (as.numeric(runtime) < 60) {
  cat("✓ Algorithm works and is reasonably fast\n")
  cat("  Full dataset will take ~", round(259659 * 6 * as.numeric(runtime) / 1000 / 3600, 1), "hours\n")
} else if (as.numeric(runtime) < 300) {
  cat("⚠️  Algorithm works but is slow\n")
  cat("  Consider optimization or parallelization\n")
} else {
  cat("❌ Algorithm is too slow for practical use\n")
  cat("  Need major optimization or different approach\n")
}
