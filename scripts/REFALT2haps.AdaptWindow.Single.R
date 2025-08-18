#!/usr/bin/env Rscript

# REFALT2haps Adaptive Window - Binary Distinguishability Check
# This script progressively expands window size until all 8 founders can be distinguished
# Outputs estimate_OK: 1 if eventually distinguishable, 0 if not (within size limits)
# 
# Usage: Rscript REFALT2haps.AdaptWindow.Single.R <chr> <parfile> <mydir> <h_cutoff>
# Example: Rscript REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 4

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript REFALT2haps.AdaptWindow.Single.R <chr> <parfile> <mydir> <h_cutoff>\n")
  cat("Example: Rscript REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 4\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2] 
mydir <- args[3]
h_cutoff_from_cmdline <- as.numeric(args[4])

# Source the parameter file
source(parfile)

# CRITICAL: Use command line h_cutoff, not parameter file value
h_cutoff <- h_cutoff_from_cmdline

# Define file paths
filein <- paste0(dirname(mydir), "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Adaptive Window - Distinguishability Check ===\n")
cat("Chromosome:", mychr, "\n")
cat("H_cutoff:", h_cutoff, "\n")
cat("Input file:", filein, "\n")
cat("Parameter file:", parfile, "\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), h_cutoff (", h_cutoff, "), samples (", length(names_in_bam), ")\n\n")

# Load data
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

cat("H_cutoff:", h_cutoff, "\n")
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
  
  # Process each sample
  for (sample_name in non_founder_samples) {
    
    # Adaptive window expansion: start small, grow until distinguishable or max size
    window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)  # Progressive sizes
    final_window_size <- NA
    estimate_OK <- 0
    final_n_snps <- 0
    
    for (window_size in window_sizes) {
      # Define window boundaries
      window_start <- test_pos - window_size/2
      window_end <- test_pos + window_size/2
      
      # Get SNPs in window for founders only (data is already quality-filtered)
      window_snps_long <- df3 %>%
        filter(POS >= window_start & POS <= window_end & name %in% founders)
      
      # Convert to wide format (rows = positions, columns = founders)
      wide_data <- window_snps_long %>%
        select(POS, name, freq) %>%
        pivot_wider(names_from = name, values_from = freq)
      
      # Check if we have all founder columns and enough data
      if (ncol(wide_data) < length(founders) + 1 || nrow(wide_data) < 10) {
        next  # Try larger window
      }
      
      # Get founder matrix (no need for additional quality filtering)
      founder_matrix <- wide_data %>%
        select(all_of(founders)) %>%
        as.matrix()
      
      # Convert to matrix for clustering
      founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]
      
      if (nrow(founder_matrix_clean) < 10) {
        next  # Try larger window
      }
      
      # Hierarchical clustering to check distinguishability
      tryCatch({
        distances <- dist(t(founder_matrix_clean), method = "euclidean")
        hclust_result <- hclust(distances, method = "ward.D2")
        
        # Cut tree at h_cutoff
        groups <- cutree(hclust_result, h = h_cutoff)
        n_groups <- length(unique(groups))
        
        # Check if all founders can be distinguished
        if (n_groups == length(founders)) {
          final_window_size <- window_size
          estimate_OK <- 1
          final_n_snps <- nrow(wide_data)
          break  # Success! Stop expanding
        }
        
      }, error = function(e) {
        # Clustering failed, try larger window
      })
    }
    
    # Record result (either successful at some window size, or failed at all sizes)
    if (is.na(final_window_size)) {
      final_window_size <- max(window_sizes)  # Used largest window but failed
      # final_n_snps should be from the last attempt
      window_start <- test_pos - final_window_size/2
      window_end <- test_pos + final_window_size/2
      window_snps_long <- df3 %>%
        filter(POS >= window_start & POS <= window_end & name %in% founders)
      wide_data <- window_snps_long %>%
        select(POS, name, freq) %>%
        pivot_wider(names_from = name, values_from = freq)
      final_n_snps <- nrow(wide_data)
    }
    
    # Now run LSEI estimation using the optimal window size
    tryCatch({
      # Get data for the final window
      window_start <- test_pos - final_window_size/2
      window_end <- test_pos + final_window_size/2
      
      window_snps_long <- df3 %>%
        filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
      
      if (nrow(window_snps_long) == 0) {
        # No data for LSEI
        na_freqs <- rep(NA, length(founders))
        names(na_freqs) <- founders
        
        result_row <- as.list(c(
          chr = mychr,
          pos = test_pos,
          sample = sample_name,
          h_cutoff = h_cutoff,
          final_window_size = final_window_size,
          n_snps = final_n_snps,
          estimate_OK = NA,
          na_freqs
        ))
      } else {
        # Prepare data for LSEI
        wide_data <- window_snps_long %>%
          select(POS, name, freq) %>%
          pivot_wider(names_from = name, values_from = freq)
        
        # Extract founder and sample data
        founder_matrix <- wide_data %>%
          select(all_of(founders)) %>%
          as.matrix()
        sample_freqs <- wide_data %>%
          pull(!!sample_name)
        
        # Remove rows with any NAs
        complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
        founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
        sample_freqs_clean <- sample_freqs[complete_rows]
        
        if (nrow(founder_matrix_clean) == 0) {
          # No complete data for LSEI
          na_freqs <- rep(NA, length(founders))
          names(na_freqs) <- founders
          
          result_row <- as.list(c(
            chr = mychr,
            pos = test_pos,
            sample = sample_name,
            h_cutoff = h_cutoff,
            final_window_size = final_window_size,
            n_snps = final_n_snps,
            estimate_OK = NA,
            na_freqs
          ))
        } else {
          # Run LSEI
          E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
          F <- 1.0
          G <- diag(length(founders))  # Non-negativity constraints
          H <- matrix(rep(0.0003, length(founders)))  # Lower bound
          
          lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                       E = E, F = F, G = G, H = H)
          
          if (lsei_result$IsError == 0) {
            # LSEI successful
            haplotype_freqs <- lsei_result$X
            names(haplotype_freqs) <- founders
            
            result_row <- as.list(c(
              chr = mychr,
              pos = test_pos,
              sample = sample_name,
              h_cutoff = h_cutoff,
              final_window_size = final_window_size,
              n_snps = final_n_snps,
              estimate_OK = estimate_OK,
              haplotype_freqs
            ))
          } else {
            # LSEI failed
            na_freqs <- rep(NA, length(founders))
            names(na_freqs) <- founders
            
            result_row <- as.list(c(
              chr = mychr,
              pos = test_pos,
              sample = sample_name,
              h_cutoff = h_cutoff,
              final_window_size = final_window_size,
              n_snps = final_n_snps,
              estimate_OK = NA,  # NA because LSEI failed
              na_freqs
            ))
          }
        }
      }
    }, error = function(e) {
      # LSEI error
      na_freqs <- rep(NA, length(founders))
      names(na_freqs) <- founders
      
      result_row <- as.list(c(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        final_window_size = final_window_size,
        n_snps = final_n_snps,
        estimate_OK = NA,  # NA because of error
        na_freqs
      ))
    })
    
    results_list[[length(results_list) + 1]] <- result_row
  }
}

# Convert results to data frame
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Save results
  output_file <- paste0(mydir, "/adaptive_window_h", h_cutoff, "_results_", mychr, ".RDS")
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
      avg_window_size = mean(final_window_size, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(distinguishable))
  print(sample_counts)
  
  # Show window size distribution for successful cases
  if (distinguishable_count > 0) {
    cat("\nWindow size distribution for successful cases:\n")
    window_summary <- results_df %>%
      filter(estimate_OK == 1) %>%
      count(final_window_size, name = "count") %>%
      mutate(percentage = count / distinguishable_count * 100) %>%
      arrange(final_window_size)
    print(window_summary)
  }
  
  cat("✓ Adaptive distinguishability analysis completed successfully\n")
} else {
  cat("\n❌ No results obtained!\n")
  quit(status = 1)
}
