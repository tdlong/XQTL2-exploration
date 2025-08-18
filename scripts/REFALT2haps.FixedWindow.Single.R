#!/usr/bin/env Rscript

# REFALT2haps Fixed Window - Binary Distinguishability Check
# This script checks if all 8 founders can be distinguished at a fixed window size
# Outputs estimate_OK: 1 if distinguishable, 0 if not
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

# Define file paths
filein <- paste0(dirname(mydir), "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Fixed Window - Distinguishability Check ===\n")
cat("Chromosome:", mychr, "\n")
cat("Window size:", window_size_kb, "kb (", window_size_bp, "bp)\n")
cat("Input file:", filein, "\n")
cat("Parameter file:", parfile, "\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), h_cutoff (", h_cutoff, "), samples (", length(names_in_bam), ")\n\n")

# Load data (same as test script)
cat("Loading data...\n")
df <- read.table(filein, header = TRUE)

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

# Transform to wide format and apply quality filter ONCE
cat("Converting to wide format and applying quality filter...\n")

# Get founder data in wide format
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

cat("Quality-filtered positions:", length(quality_filtered_positions), "\n")

# Keep only high-quality positions in the full dataset
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

# Use samples defined in parameter file instead of auto-detecting
non_founder_samples <- names_in_bam

cat("Using samples from parameter file:", length(non_founder_samples), "\n")
cat("Samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Define scanning positions (500kb to end-500kb, 10kb steps)
chromosome_length <- max(df$POS, na.rm = TRUE)
if (!is.finite(chromosome_length)) {
  chromosome_length <- max(euchromatin_2R)  # fallback to euchromatin
}

scan_start <- 500000
scan_end <- chromosome_length - 500000
scan_positions <- seq(scan_start, scan_end, by = step)
cat("Using step size:", step, "bp\n")

cat("Chromosome length:", chromosome_length, "bp\n")
cat("Scanning from:", scan_start, "to", scan_end, "bp\n")
cat("Total positions to scan:", length(scan_positions), "\n\n")

# Initialize results table
results_list <- list()

# Scan each position
for (pos_idx in seq_along(scan_positions)) {
  test_pos <- scan_positions[pos_idx]
  
  if (pos_idx %% 100 == 0) {
    cat("Processing position", pos_idx, "of", length(scan_positions), "(", test_pos, "bp)\n")
  }
  
  # Define window boundaries
  window_start <- test_pos - window_size_bp/2
  window_end <- test_pos + window_size_bp/2
  
  # Process each sample
  for (sample_name in non_founder_samples) {
    
    # Get SNPs in window for founders only (data is already quality-filtered)
    window_snps_long <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    # Convert to wide format (rows = positions, columns = founders + sample)
    wide_data <- window_snps_long %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    # Check if we have all required columns and enough data
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) {
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        window_size = window_size_bp,
        n_snps = nrow(wide_data),
        estimate_OK = 0  # Cannot distinguish - missing founders or insufficient data
      )
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    # Get founder matrix (no need for additional quality filtering)
    founder_matrix <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    
    # Convert to matrix for clustering
    founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]
    
    if (nrow(founder_matrix_clean) < 10) {
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        window_size = window_size_bp,
        n_snps = nrow(wide_data),
        estimate_OK = 0  # Cannot distinguish - insufficient clean data
      )
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    # Hierarchical clustering to check distinguishability
    tryCatch({
      distances <- dist(t(founder_matrix_clean), method = "euclidean")
      hclust_result <- hclust(distances, method = "ward.D2")
      
      # Run LSEI first to get actual haplotype frequency estimates
      tryCatch({
        # LSEI constraints: sum to 1, non-negative
        E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
        F <- 1.0
        G <- diag(length(founders))  # Non-negativity constraints
        H <- matrix(rep(0.0003, length(founders)))  # Lower bound
        
        lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs, 
                                     E = E, F = F, G = G, H = H)
        
        if (lsei_result$IsError == 0) {
          # LSEI successful - now check distinguishability
          haplotype_freqs <- lsei_result$X
          names(haplotype_freqs) <- founders
          
          # Check distinguishability using clustering
          groups <- cutree(hclust_result, h = h_cutoff)
          n_groups <- length(unique(groups))
          estimate_OK <- ifelse(n_groups == length(founders), 1, 0)
          
          result_row <- as.list(c(
            chr = mychr,
            pos = test_pos,
            sample = sample_name,
            window_size = window_size_bp,
            n_snps = nrow(wide_data),
            estimate_OK = estimate_OK,
            haplotype_freqs
          ))
        } else {
          # LSEI failed - output NAs for frequencies and estimate_OK
          na_freqs <- rep(NA, length(founders))
          names(na_freqs) <- founders
          
          result_row <- as.list(c(
            chr = mychr,
            pos = test_pos,
            sample = sample_name,
            window_size = window_size_bp,
            n_snps = nrow(wide_data),
            estimate_OK = NA,  # NA because LSEI failed
            na_freqs
          ))
        }
      }, error = function(e) {
        # LSEI error - output NAs for frequencies
        na_freqs <- rep(NA, length(founders))
        names(na_freqs) <- founders
        
        result_row <- as.list(c(
          chr = mychr,
          pos = test_pos,
          sample = sample_name,
          window_size = window_size_bp,
          n_snps = nrow(wide_data),
          estimate_OK = NA,  # NA because LSEI failed
          na_freqs
        ))
      })
      
      results_list[[length(results_list) + 1]] <- result_row
      
    }, error = function(e) {
      # Cannot distinguish - clustering failed, output NAs for all frequencies
      na_freqs <- rep(NA, length(founders))
      names(na_freqs) <- founders
      
      result_row <- as.list(c(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        window_size = window_size_bp,
        n_snps = nrow(wide_data),
        estimate_OK = NA,  # NA because clustering failed
        na_freqs
      ))
      
      results_list[[length(results_list) + 1]] <- result_row
    })
  }
}

# Convert results to data frame
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Save results
  output_file <- paste0(mydir, "/fixed_window_", window_size_kb, "kb_results_", mychr, ".RDS")
  saveRDS(results_df, output_file)
  
  cat("\n=== RESULTS SUMMARY ===\n")
  cat("Total position/sample combinations:", nrow(results_df), "\n")
  cat("Output file:", output_file, "\n")
  cat("File size:", file.size(output_file), "bytes\n")
  
  # Calculate distinguishability rate
  total_combinations <- length(scan_positions) * length(non_founder_samples)
  distinguishable_count <- sum(results_df$estimate_OK)
  distinguishability_rate <- distinguishable_count / total_combinations * 100
  
  cat("\nDistinguishability rate:", round(distinguishability_rate, 1), "% (", distinguishable_count, "of", total_combinations, "combinations)\n")
  cat("Positions scanned:", length(scan_positions), "\n")
  cat("Samples processed:", length(non_founder_samples), "\n")
  
  # Show sample summary
  cat("\nDistinguishability per sample:\n")
  sample_counts <- results_df %>%
    group_by(sample) %>%
    summarize(
      total_results = n(),
      distinguishable = sum(estimate_OK),
      distinguishability_rate = distinguishable / length(scan_positions) * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(distinguishable))
  print(sample_counts)
  
  cat("✓ Distinguishability analysis completed successfully\n")
} else {
  cat("\n❌ No results obtained!\n")
  quit(status = 1)
}
