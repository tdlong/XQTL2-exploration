#!/usr/bin/env Rscript

# =============================================================================
# REFALT2haps Parameter Assessment Script - KISS Version
# =============================================================================
# Simple script to test how window size affects haplotype frequency estimates
# Usage: Rscript scripts/REFALT2haps.Assessment.R chr parfile mydir start_pos end_pos

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript REFALT2haps.Assessment.R chr parfile mydir start_pos end_pos\n")
  cat("Example: Rscript REFALT2haps.Assessment.R chr2L helpfiles/haplotype_parameters.R process/test 18000000 20000000\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
start_pos <- as.numeric(args[4])
end_pos <- as.numeric(args[5])

# Calculate midpoint
midpoint <- (start_pos + end_pos) / 2

# Source the parameter file
source(parfile)

# Define file paths
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Fixed Window Assessment (KISS) ===\n")
cat("Chromosome:", mychr, "\n")
cat("Region:", start_pos, "-", end_pos, "bp\n")
cat("Midpoint:", midpoint, "bp\n")
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

# Define window sizes to test
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000, 1000000)
results <- matrix(NA, nrow = length(founders), ncol = length(window_sizes))
rownames(results) <- founders
colnames(results) <- window_sizes

cat("Testing window sizes:", paste(window_sizes/1000, "kb", collapse = ", "), "\n")
cat("Using fixed h_cutoff:", h_cutoff, "\n\n")

# Test first sample
test_sample <- names_in_bam[1]

for (i in seq_along(window_sizes)) {
  ws <- window_sizes[i]
  
  # Create window around midpoint
  window_start <- max(0, midpoint - ws/2)
  window_end <- midpoint + ws/2
  
  cat("\n--- Testing window size:", ws/1000, "kb ---\n")
  cat("Window:", window_start, "-", window_end, "bp\n")
  
  # Filter SNPs in window and reshape data
  window_snps <- df3 %>%
    filter(CHROM == mychr &
           POS >= window_start & 
           POS <= window_end &
           (name %in% founders | name == test_sample)) %>%
    select(-c(CHROM, N)) %>%
    pivot_wider(names_from = name, values_from = freq) %>%
    filter(!is.na(!!sym(test_sample)))
  
  cat("SNPs in window:", nrow(window_snps), "\n")
  
  if (nrow(window_snps) > 0) {
    # Extract founder matrix and sample frequencies
    founder_matrix <- window_snps %>% select(all_of(founders))
    sample_freqs <- window_snps[[test_sample]]
    
    # Convert to matrix for clustering
    founder_matrix <- as.matrix(founder_matrix)
    
    # Cluster founders based on similarity
    founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
    
    # Count groups and show composition
    n_groups <- length(unique(founder_clusters))
    cat("Founder groups:", n_groups, "\n")
    
    for (group_id in unique(founder_clusters)) {
      group_founders <- names(founder_clusters[founder_clusters == group_id])
      cat("  Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
    }
    
    # Build constraint matrix
    n_founders <- ncol(founder_matrix)
    E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    
    # Add group constraints for each cluster (only if multiple groups)
    unique_clusters <- unique(founder_clusters)
    if (n_groups > 1) {
      for (cluster_id in unique_clusters) {
        cluster_founders <- which(founder_clusters == cluster_id)
        if (length(cluster_founders) > 1) {
          constraint_row <- rep(0, n_founders)
          constraint_row[cluster_founders] <- 1
          E <- rbind(E, constraint_row)
          F <- c(F, 1.0)
        }
      }
    }
    
    cat("Constraints: Sum to 1 +", nrow(E)-1, "group constraints\n")
    
    # Solve constrained least squares with proper bounds
    tryCatch({
      # G = diag(n_founders) sets lower bounds (>= 0.0003)
      # H = matrix(rep(0.0003, n_founders)) sets minimum frequency
      result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F,
                    G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
      
      cat("✓ Estimation successful\n")
      cat("Frequency estimates:\n")
      for (j in seq_along(founders)) {
        cat("  ", founders[j], ":", sprintf("%.4f", result$X[j]), "\n")
      }
      
      results[, i] <- result$X
    }, error = function(e) {
      cat("✗ Estimation failed:", e$message, "\n")
      results[, i] <- rep(NA, length(founders))
    })
  } else {
    cat("No SNPs found in window\n")
  }
}

# Output results table
cat("Haplotype Frequency Estimates vs Window Size\n")
cat("Sample:", test_sample, "\n")
cat("Position:", midpoint, "bp\n")
cat("Fixed h_cutoff:", h_cutoff, "\n\n")

# Print table
cat("Founder\t", paste(window_sizes/1000, "kb", collapse = "\t"), "\n", sep = "")
for (i in seq_along(founders)) {
  cat(founders[i], "\t", paste(sprintf("%.4f", results[i, ]), collapse = "\t"), "\n", sep = "")
}

cat("\n=== Done ===\n")
