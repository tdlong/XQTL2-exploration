#!/usr/bin/env Rscript

# =============================================================================
# REFALT2haps Adaptive Window Testing Script - KISS Version
# =============================================================================
# Simple script to test how h_cutoff affects haplotype frequency estimates
# Usage: Rscript scripts/REFALT2haps.AdaptWindow.R chr parfile mydir test_pos test_window_size

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript REFALT2haps.AdaptWindow.R chr parfile mydir test_pos test_window_size\n")
  cat("Example: Rscript REFALT2haps.AdaptWindow.R chr3R helpfiles/haplotype_parameters.R process/test 10000000 2000000\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
test_pos <- as.numeric(args[4])
test_window_size <- as.numeric(args[5])

# Source the parameter file
source(parfile)

# Define file paths
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Adaptive Window Test (KISS) ===\n")
cat("Chromosome:", mychr, "\n")
cat("Test position:", test_pos, "\n")
cat("Input file:", filein, "\n\n")

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

# Subset dataset to only high-quality SNPs
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Create window around test position (use a reasonable window size, not the full test_window_size)
base_window_size <- 1000000  # 1Mb base window
window_start <- max(0, test_pos - base_window_size/2)
window_end <- test_pos + base_window_size/2

# Filter SNPs in window
window_snps <- df3 %>%
  filter(CHROM == mychr &
         POS > window_start &
         POS < window_end &
         (name %in% founders | name %in% names_in_bam)) %>%
  select(-c(CHROM, N)) %>%
  pivot_wider(names_from = name, values_from = freq) %>%
  pivot_longer(!c("POS", matches(founders)), 
              names_to = "sample", values_to = "freq") %>%
  select(-POS)

if (nrow(window_snps) == 0) {
  stop("No SNPs found in window!")
}

# Test different h_cutoffs (simplified for debugging)
h_cutoffs <- c(2, 10, 100)  # Low, middle, high
window_sizes <- c(10000, 25000, 50000, 100000, 200000, 500000, 1000000)

cat("Testing h_cutoffs:", paste(h_cutoffs, collapse = ", "), "\n")
cat("Window sizes:", paste(window_sizes/1000, "kb", collapse = " → "), "\n")
cat("Test position:", test_pos, "bp\n\n")

# Test first sample
test_sample <- names_in_bam[1]

# Store final results for each h_cutoff
final_results <- matrix(NA, nrow = length(founders), ncol = length(h_cutoffs))
rownames(final_results) <- founders
colnames(final_results) <- h_cutoffs

# Test each h_cutoff
for (hc_idx in seq_along(h_cutoffs)) {
  hc <- h_cutoffs[hc_idx]
  cat("\n=== Testing h_cutoff:", hc, "===\n")
  
  previous_constraints <- NULL
  
  # Run adaptive window algorithm for this h_cutoff
  for (window_idx in seq_along(window_sizes)) {
    current_window_size <- window_sizes[window_idx]
    cat("\n--- Window size:", current_window_size/1000, "kb ---\n")
    
    # Create window around test position (this is the adaptive part)
    window_start <- max(0, test_pos - current_window_size/2)
    window_end <- test_pos + current_window_size/2
    cat("Window:", window_start, "-", window_end, "bp (centered at", test_pos, "bp)\n")
    
    # Filter SNPs in this expanding window
    window_snps_current <- df3 %>%
      filter(CHROM == mychr &
             POS > window_start &
             POS < window_end &
             (name %in% founders | name %in% names_in_bam)) %>%
      select(-c(CHROM, N)) %>%
      pivot_wider(names_from = name, values_from = freq) %>%
      pivot_longer(!c("POS", matches(founders)), 
                  names_to = "sample", values_to = "freq") %>%
      select(-POS)
    
    # Count SNPs in this window
    n_snps <- df3 %>%
      filter(CHROM == mychr &
             POS > window_start &
             POS < window_end &
             name %in% founders) %>%
      distinct(POS) %>%
      nrow()
    
    cat("SNPs in window:", n_snps, "\n")
    
    if (nrow(window_snps_current) == 0) {
      cat("No SNPs found in window\n")
      next
    }
    
    # Get sample data for current window
    sample_data_current <- window_snps_current %>% filter(sample == test_sample)
    
    # Extract founder matrix and sample frequencies
    founder_matrix <- sample_data_current %>% select(matches(founders))
    sample_freqs <- sample_data_current$freq
    
    # Filter for non-NA values
    valid_positions <- !is.na(sample_freqs)
    sample_freqs <- sample_freqs[valid_positions]
    founder_matrix <- founder_matrix[valid_positions, ]
    
    if (nrow(founder_matrix) > 0) {
      # Convert to matrix for clustering
      founder_matrix <- as.matrix(founder_matrix)
      
      # Cluster founders based on similarity
      founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = hc)
      
      # Count groups
      n_groups <- length(unique(founder_clusters))
      cat("Founder groups:", n_groups, "\n")
      
      # Show group composition
      for (group_id in unique(founder_clusters)) {
        group_founders <- names(founder_clusters[founder_clusters == group_id])
        cat("  Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
      }
      
      # Build constraint matrix
      n_founders <- ncol(founder_matrix)
      E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      
      # Add group constraints for each cluster (only for groups with multiple founders)
      unique_clusters <- unique(founder_clusters)
      multi_founder_groups <- 0
      for (cluster_id in unique_clusters) {
        cluster_founders <- which(founder_clusters == cluster_id)
        if (length(cluster_founders) > 1) {
          constraint_row <- rep(0, n_founders)
          constraint_row[cluster_founders] <- 1
          E <- rbind(E, constraint_row)
          F <- c(F, 1.0)
          multi_founder_groups <- multi_founder_groups + 1
        }
      }
      
      # Solve constrained least squares with proper bounds
      tryCatch({
        # G = diag(n_founders) sets lower bounds (>= 0.0003)
        # H = matrix(rep(0.0003, n_founders)) sets minimum frequency
        result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                      G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
        
        # Show frequency estimates for this window
        cat("Frequency estimates:\n")
        for (i in seq_along(founders)) {
          cat("  ", founders[i], ":", sprintf("%.4f", result$X[i]), "\n")
        }
        
        # Store result for this h_cutoff (store when we succeed)
        final_results[, hc_idx] <- result$X
        
        # Check if all founders are separated
        if (n_groups == length(founders)) {
          cat("✓ All founders separated! Stopping for this h_cutoff.\n")
          break
        }
        
      }, error = function(e) {
        cat("Error in estimation\n")
      })
    } else {
      cat("No valid SNPs in window\n")
    }
  }
  
  cat("=== Completed h_cutoff", hc, "===\n")
}

# Output final results table
cat("\n\n=== FINAL RESULTS TABLE ===\n")
cat("Haplotype Frequency Estimates vs H_cutoff\n")
cat("Sample:", test_sample, "\n")
cat("Window:", window_start, "-", window_end, "bp\n\n")

# Print table
cat("Founder\t", paste(h_cutoffs, collapse = "\t"), "\n", sep = "")
for (i in seq_along(founders)) {
  cat(founders[i], "\t", paste(sprintf("%.4f", final_results[i, ]), collapse = "\t"), "\n", sep = "")
}

cat("\n=== Done ===\n")
