#!/usr/bin/env Rscript

# =============================================================================
# Haplotype Estimator Evaluation Script
# =============================================================================
# Evaluates 12 different haplotype estimation methods by comparing
# observed SNP frequencies to imputed frequencies from haplotype estimates
# Usage: Rscript scripts/evaluate_haplotype_estimators.R chr parfile mydir

library(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript evaluate_haplotype_estimators.R chr parfile mydir\n")
  cat("Example: Rscript evaluate_haplotype_estimators.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]

# Source the parameter file
source(parfile)

cat("=== Haplotype Estimator Evaluation ===\n")
cat("Chromosome:", mychr, "\n")
cat("Parameter file:", parfile, "\n")
cat("Output directory:", mydir, "\n\n")

# Load haplotype estimation results
cat("Loading haplotype estimation results...\n")
fixed_results_file <- paste0(mydir, "/fixed_window_results_", mychr, ".RDS")
adaptive_results_file <- paste0(mydir, "/adaptive_window_results_", mychr, ".RDS")

if (!file.exists(fixed_results_file)) {
  stop("Fixed window results file not found: ", fixed_results_file)
}
if (!file.exists(adaptive_results_file)) {
  stop("Adaptive window results file not found: ", adaptive_results_file)
}

fixed_results <- readRDS(fixed_results_file)
adaptive_results <- readRDS(adaptive_results_file)

cat("✓ Fixed window results loaded:", nrow(fixed_results), "estimates\n")
cat("✓ Adaptive window results loaded:", nrow(adaptive_results), "estimates\n")

# Load observed SNP data from REFALT
cat("Loading observed SNP data...\n")
refalt_file <- paste0(mydir, "/RefAlt.", mychr, ".txt")
if (!file.exists(refalt_file)) {
  stop("REFALT file not found: ", refalt_file)
}

df <- read.table(refalt_file, header = TRUE)

# Transform REF/ALT counts to frequencies
cat("Converting counts to frequencies...\n")
cat("Columns in REFALT data:", paste(names(df), collapse = ", "), "\n")
cat("Number of columns:", ncol(df), "\n")

# Check if we have the expected structure
if (ncol(df) < 3) {
  stop("REFALT file has insufficient columns. Expected: CHROM, POS, and sample columns")
}

# Get sample columns (exclude CHROM and POS)
sample_cols <- names(df)[!names(df) %in% c("CHROM", "POS")]
cat("Sample columns:", paste(sample_cols, collapse = ", "), "\n")

df2 <- df %>%
  pivot_longer(all_of(sample_cols), names_to = "lab", values_to = "count") %>%
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

# Filter for high-quality SNPs (same criteria as estimation scripts)
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

# Get valid SNPs for evaluation
valid_snps <- good_snps %>%
  left_join(df2, multiple = "all")

cat("✓ Valid SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n")
cat("Columns in valid_snps:", paste(names(valid_snps), collapse = ", "), "\n")
cat("Sample of valid_snps data:\n")
print(head(valid_snps, 3))

# Function to interpolate haplotype frequencies to SNP positions
interpolate_haplotype_frequencies <- function(haplotype_results, snp_positions, founders) {
  # Debug: Check the structure of haplotype_results
  cat("Debug: haplotype_results columns:", paste(names(haplotype_results), collapse = ", "), "\n")
  cat("Debug: haplotype_results structure:\n")
  print(str(haplotype_results))
  cat("Debug: Expected founders:", paste(founders, collapse = ", "), "\n")
  
  # Get unique positions where we have haplotype estimates
  haplotype_positions <- unique(haplotype_results$pos)
  haplotype_positions <- sort(haplotype_positions)
  
  # Initialize results
  interpolated_results <- list()
  
  for (snp_pos in snp_positions) {
    # Find nearest haplotype estimates (left and right)
    left_pos <- max(haplotype_positions[haplotype_positions <= snp_pos])
    right_pos <- min(haplotype_positions[haplotype_positions >= snp_pos])
    
    if (left_pos == right_pos) {
      # SNP is exactly at a haplotype estimate position
      haplotype_freqs <- haplotype_results %>%
        filter(pos == snp_pos) %>%
        select(founder, freq) %>%
        pivot_wider(names_from = founder, values_from = freq)
      
      # Debug: Check what columns we actually have after pivot_wider
      cat("Debug: After pivot_wider, columns:", paste(names(haplotype_freqs), collapse = ", "), "\n")
      
      # Only select founder columns that actually exist
      existing_founders <- intersect(founders, names(haplotype_freqs))
      if (length(existing_founders) > 0) {
        haplotype_freqs <- haplotype_freqs %>% select(all_of(existing_founders))
        interpolated_results[[as.character(snp_pos)]] <- haplotype_freqs
      } else {
        cat("Warning: No founder columns found for position", snp_pos, "\n")
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
        
        # Interpolation weight
        alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
        
        # Linear interpolation: H{snp} = α*H{left} + (1-α)*H{right}
        interpolated_freqs <- alpha * left_freqs + (1 - alpha) * right_freqs
        interpolated_results[[as.character(snp_pos)]] <- interpolated_freqs
      } else {
        cat("Warning: No founder columns found for interpolation at position", snp_pos, "\n")
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

# Function to create evaluation table
create_evaluation_table <- function(valid_snps, fixed_results, adaptive_results, founders) {
  cat("Creating evaluation table...\n")
  
  # Get unique SNP positions
  snp_positions <- valid_snps %>% distinct(CHROM, POS) %>% pull(POS)
  
  # Get unique samples (non-founders)
  all_samples <- unique(valid_snps$name)
  non_founder_samples <- all_samples[!all_samples %in% founders]
  
  # Initialize results table
  evaluation_table <- list()
  
  # Process each sample separately
  for (sample_name in non_founder_samples) {
    cat("Processing sample:", sample_name, "\n")
    
    # Filter results for this sample
    sample_fixed_results <- fixed_results %>% filter(sample == sample_name)
    sample_adaptive_results <- adaptive_results %>% filter(sample == sample_name)
    
    # Interpolate haplotype frequencies for this sample
    fixed_interpolated <- interpolate_haplotype_frequencies(sample_fixed_results, snp_positions, founders)
    adaptive_interpolated <- interpolate_haplotype_frequencies(sample_adaptive_results, snp_positions, founders)
    
    # Process each SNP position
    for (snp_pos in snp_positions) {
      # Get founder states for this SNP
      snp_founder_data <- valid_snps %>%
        filter(POS == snp_pos, name %in% founders) %>%
        arrange(match(name, founders))  # Ensure founders are in correct order
      
      if (nrow(snp_founder_data) == length(founders)) {
        # Get observed frequency for this sample at this SNP
        observed_data <- valid_snps %>%
          filter(POS == snp_pos, name == sample_name)
        
        if (nrow(observed_data) > 0) {
          observed_freq <- observed_data$freq[1]
          
          # Get founder states (0 or 1)
          founder_states <- snp_founder_data$freq
          
          # Calculate imputed frequencies for each estimator
          snp_results <- list(
            chr = mychr,
            pos = snp_pos,
            sample = sample_name,
            observed_freq = observed_freq
          )
          
          # Fixed window estimators
          window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
          for (ws in window_sizes) {
            if (as.character(snp_pos) %in% names(fixed_interpolated)) {
              haplotype_freqs <- as.numeric(fixed_interpolated[[as.character(snp_pos)]])
              imputed_freq <- calculate_imputed_snp_frequency(haplotype_freqs, founder_states)
              snp_results[[paste0("fixed_", ws/1000, "kb")]] <- imputed_freq
            } else {
              snp_results[[paste0("fixed_", ws/1000, "kb")]] <- NA
            }
          }
          
          # Adaptive window estimators
          h_cutoffs <- c(2, 4, 6, 8, 10, 20)
          for (hc in h_cutoffs) {
            if (as.character(snp_pos) %in% names(adaptive_interpolated)) {
              haplotype_freqs <- as.numeric(adaptive_interpolated[[as.character(snp_pos)]])
              imputed_freq <- calculate_imputed_snp_frequency(haplotype_freqs, founder_states)
              snp_results[[paste0("adaptive_h", hc)]] <- imputed_freq
            } else {
              snp_results[[paste0("adaptive_h", hc)]] <- NA
            }
          }
          
          evaluation_table[[length(evaluation_table) + 1]] <- snp_results
        }
      }
    }
  }
  
  # Convert to data frame
  if (length(evaluation_table) > 0) {
    evaluation_df <- bind_rows(evaluation_table)
    return(evaluation_df)
  } else {
    return(NULL)
  }
}

# Function to calculate MSE for each estimator
calculate_mse_by_estimator <- function(evaluation_table) {
  cat("Calculating MSE for each estimator...\n")
  
  # Get estimator column names (exclude chr, pos, sample, observed_freq)
  estimator_cols <- names(evaluation_table)[!names(evaluation_table) %in% c("chr", "pos", "sample", "observed_freq")]
  
  mse_results <- list()
  
  for (estimator in estimator_cols) {
    # Calculate MSE: mean((observed - imputed)^2)
    mse <- evaluation_table %>%
      filter(!is.na(!!sym(estimator))) %>%
      summarize(
        estimator = estimator,
        mse = mean((observed_freq - !!sym(estimator))^2, na.rm = TRUE),
        n_comparisons = n(),
        mean_observed = mean(observed_freq, na.rm = TRUE),
        mean_imputed = mean(!!sym(estimator), na.rm = TRUE)
      )
    
    mse_results[[length(mse_results) + 1]] <- mse
  }
  
  mse_df <- bind_rows(mse_results)
  return(mse_df)
}

# Main execution
cat("Starting evaluation...\n")

# Create evaluation table
evaluation_table <- create_evaluation_table(valid_snps, fixed_results, adaptive_results, founders)

if (!is.null(evaluation_table)) {
  cat("✓ Evaluation table created:", nrow(evaluation_table), "rows\n")
  
  # Calculate MSE
  mse_results <- calculate_mse_by_estimator(evaluation_table)
  
  # Save results
  evaluation_file <- paste0(mydir, "/evaluation_table_", mychr, ".RDS")
  mse_file <- paste0(mydir, "/mse_results_", mychr, ".RDS")
  
  saveRDS(evaluation_table, evaluation_file)
  saveRDS(mse_results, mse_file)
  
  cat("✓ Results saved:\n")
  cat("  Evaluation table:", evaluation_file, "\n")
  cat("  MSE results:", mse_file, "\n")
  
  # Display MSE results
  cat("\n=== MSE Results (Ranked by Performance) ===\n")
  mse_ranked <- mse_results %>% arrange(mse)
  print(mse_ranked)
  
  cat("\n=== Evaluation Complete ===\n")
  cat("Best estimator:", mse_ranked$estimator[1], "(MSE =", sprintf("%.6f", mse_ranked$mse[1]), ")\n")
  
} else {
  cat("✗ No evaluation data could be created\n")
}
