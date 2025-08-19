#!/usr/bin/env Rscript

# REFALT2haps Fixed Window - Haplotype Frequency Estimation
# This script estimates haplotype frequencies using a fixed window size
# and outputs both frequency estimates AND distinguishability quality flag
# 
# Usage: Rscript REFALT2haps.FixedWindow.Single.R <chr> <parfile> <mydir> <window_size_kb>
# Example: Rscript REFALT2haps.FixedWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 50

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript REFALT2haps.FixedWindow.Single.R <chr> <parfile> <mydir> <window_size_kb>\n")
  cat("Example: Rscript REFALT2haps.FixedWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 50\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
window_size_kb <- as.numeric(args[4])

# Convert kb to bp
window_size_bp <- window_size_kb * 1000

# Source the parameter file
source(parfile)

# Display parameters loaded from parameter file
cat("Parameters loaded from", parfile, ":\n")
cat("  Founders (", length(founders), "):", paste(founders, collapse=", "), "\n")
cat("  Step:", step, "\n")
cat("  Samples to process (", length(names_in_bam), "):", paste(names_in_bam, collapse=", "), "\n")

# Load and transform data (same as test script)
filein <- file.path(mydir, paste0("RefAlt.", mychr, ".txt"))
df <- read.table(filein, header = TRUE)

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

# Transform to wide format and apply quality filter ONCE
founder_wide <- df2 %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

# Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
quality_filtered_positions <- founder_wide %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  ) %>%
  pull(POS)

# Filter to quality positions and include sample data
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

# Get non-founder samples from parameter file
non_founder_samples <- names_in_bam

# Set up scanning positions based on step parameter
scan_positions <- seq(from = min(quality_filtered_positions), 
                     to = max(quality_filtered_positions), 
                     by = step)

# Initialize results list
results_list <- list()

# Process each position and sample combination
for (test_pos in scan_positions) {
  for (sample_name in non_founder_samples) {
    
    # Calculate window boundaries
    window_start <- test_pos - window_size_bp/2
    window_end <- test_pos + window_size_bp/2
    
    # Get SNPs in window for both founders AND sample (data is already quality-filtered)
    window_snps_long <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    # Skip if no data
    if (nrow(window_snps_long) == 0) {
      next
    }
    
    # Convert to wide format (rows = positions, columns = founders + sample)
    wide_data <- window_snps_long %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    # Check if we have all founder columns + sample and enough data
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) {
      next
    }
    
    # Get founder matrix and sample frequencies for LSEI
    founder_matrix <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    
    sample_freqs <- wide_data %>%
      pull(!!sample_name)
    
    # Remove rows with any NAs (both founder and sample data must be complete)
    complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    
    # Skip if insufficient clean data
    if (nrow(founder_matrix_clean) < 10) {
      next
    }
    
    # Initialize result variables
    estimate_OK <- NA
    haplotype_freqs <- rep(NA, length(founders))
    names(haplotype_freqs) <- founders
    
    # Run LSEI to get actual haplotype frequency estimates
    tryCatch({
      # LSEI constraints: sum to 1, non-negative
      E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      G <- diag(length(founders))  # Non-negativity constraints
      H <- matrix(rep(0.0003, length(founders)))  # Lower bound
      
      lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                   E = E, F = F, G = G, H = H)
      
      if (lsei_result$IsError == 0) {
        # LSEI successful - get frequencies
        haplotype_freqs <- lsei_result$X
        names(haplotype_freqs) <- founders
        
        # Check distinguishability using clustering
        distances <- dist(t(founder_matrix_clean), method = "euclidean")
        hclust_result <- hclust(distances, method = "ward.D2")
        groups <- cutree(hclust_result, h = h_cutoff)
        n_groups <- length(unique(groups))
        
        # Set estimate_OK based on distinguishability
        estimate_OK <- ifelse(n_groups == length(founders), 1, 0)
        
      } else {
        # LSEI failed
        estimate_OK <- NA
      }
    }, error = function(e) {
      # LSEI error
      estimate_OK <- NA
    })
    
    # CREATE result_row STRUCTURE (same as test script)
    result_row <- list(
      chr = mychr,
      pos = test_pos,
      sample = sample_name,
      window_size = window_size_bp,
      n_snps = nrow(wide_data),
      estimate_OK = estimate_OK
    )
    
    # Add founder frequencies as named columns
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- haplotype_freqs[i]
    }
    
    results_list[[length(results_list) + 1]] <- result_row
  }
}

# Convert to data frame and save
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(mydir, "haplotype_results")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  output_file <- file.path(output_dir, paste0("fixed_window_", window_size_kb, "kb_results_", mychr, ".RDS"))
  saveRDS(results_df, output_file)
  
  cat("Fixed window haplotype estimation complete.\n")
  cat("Processed", nrow(results_df), "position/sample combinations\n")
  cat("Results saved to:", output_file, "\n")
  
} else {
  cat("No results generated - check data and parameters\n")
  quit(status = 1)
}