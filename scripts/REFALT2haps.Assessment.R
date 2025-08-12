#!/usr/bin/env Rscript

# =============================================================================
# REFALT2haps Parameter Assessment Script - Chromosome Scanner
# =============================================================================
# Script to assess how window size affects haplotype estimation across entire chromosome
# Usage: Rscript scripts/REFALT2haps.Assessment.R chr parfile mydir [verbose]

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript REFALT2haps.Assessment.R chr parfile mydir [verbose]\n")
  cat("Example: Rscript REFALT2haps.Assessment.R chr2L helpfiles/haplotype_parameters.R process/test\n")
  cat("Add 'verbose' as 4th argument for detailed output\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
verbose <- ifelse(length(args) >= 4 && args[4] == "verbose", TRUE, FALSE)

# Source the parameter file
source(parfile)

# Define file paths
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Fixed Window Assessment - Chromosome Scanner ===\n")
cat("Chromosome:", mychr, "\n")
cat("Input file:", filein, "\n")
cat("Verbose mode:", verbose, "\n\n")

# Load data
if (verbose) cat("Loading data...\n")
df <- read.table(filein, header = TRUE)

# Transform REF/ALT counts to frequencies
if (verbose) cat("Converting counts to frequencies...\n")
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
if (verbose) cat("Filtering for high-quality SNPs...\n")
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

# Define window sizes to test (removed 1000kb as requested)
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)

# Get all non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

if (verbose) {
  cat("Testing window sizes:", paste(window_sizes/1000, "kb", collapse = ", "), "\n")
  cat("Using fixed h_cutoff:", h_cutoff, "\n")
  cat("Non-founder samples:", length(non_founder_samples), "\n")
  cat("Samples:", paste(non_founder_samples, collapse = ", "), "\n\n")
}

# Define scanning positions (500kb to end-500kb, 10kb steps)
chromosome_length <- max(df$POS)
scan_start <- 500000
scan_end <- chromosome_length - 500000
scan_positions <- seq(scan_start, scan_end, by = 10000)

if (verbose) {
  cat("Chromosome length:", chromosome_length, "bp\n")
  cat("Scanning from:", scan_start, "to", scan_end, "bp\n")
  cat("Total positions to scan:", length(scan_positions), "\n\n")
}

# Initialize results table
results_list <- list()

# Scan each position
for (pos_idx in seq_along(scan_positions)) {
  test_pos <- scan_positions[pos_idx]
  
  if (verbose) {
    cat("Position", pos_idx, "/", length(scan_positions), ":", test_pos, "bp\n")
  }
  
  # Test each window size for this position
  for (ws_idx in seq_along(window_sizes)) {
    ws <- window_sizes[ws_idx]
    
    # Create window around test position
    window_start <- max(0, test_pos - ws/2)
    window_end <- test_pos + ws/2
    
    # Filter SNPs in this window
    window_snps <- df3 %>%
      filter(CHROM == mychr &
             POS > window_start &
             POS < window_end &
             (name %in% founders | name %in% non_founder_samples)) %>%
      select(-c(CHROM, N)) %>%
      pivot_wider(names_from = name, values_from = freq) %>%
      pivot_longer(!c("POS", matches(founders)), 
                  names_to = "sample", values_to = "freq") %>%
      select(-POS)
    
    if (nrow(window_snps) > 0) {
      # Test each non-founder sample
      for (sample_name in non_founder_samples) {
        sample_data <- window_snps %>% filter(sample == sample_name)
        
        if (nrow(sample_data) > 0) {
          # Extract founder matrix and sample frequencies
          founder_matrix <- sample_data %>% select(matches(founders))
          sample_freqs <- sample_data$freq
          
          # Filter for non-NA values
          valid_positions <- !is.na(sample_freqs)
          sample_freqs <- sample_freqs[valid_positions]
          founder_matrix <- founder_matrix[valid_positions, ]
          
          if (nrow(founder_matrix) > 0) {
            # Convert to matrix for clustering
            founder_matrix <- as.matrix(founder_matrix)
            
            # Cluster founders based on similarity
            founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
            
            # Count groups
            n_groups <- length(unique(founder_clusters))
            
            # Build constraint matrix (simple: sum to 1 + individual bounds)
            n_founders <- ncol(founder_matrix)
            E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
            F <- 1.0
            
            # Solve constrained least squares
            tryCatch({
              result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                            G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
              
              # Store results
              for (founder_idx in seq_along(founders)) {
                founder_name <- founders[founder_idx]
                freq_estimate <- result$X[founder_idx]
                
                results_list[[length(results_list) + 1]] <- list(
                  chr = mychr,
                  pos = test_pos,
                  sample = sample_name,
                  window_size = ws,
                  founder = founder_name,
                  freq = freq_estimate
                )
              }
              
              if (verbose) {
                cat("  ✓", sample_name, "window", ws/1000, "kb:", n_groups, "groups, estimation successful\n")
              }
              
            }, error = function(e) {
              if (verbose) {
                cat("  ✗", sample_name, "window", ws/1000, "kb:", n_groups, "groups, estimation failed:", e$message, "\n")
              }
            })
          }
        }
      }
    }
  }
  
  # Progress indicator
  if (pos_idx %% 100 == 0) {
    cat("Progress:", pos_idx, "/", length(scan_positions), "positions completed\n")
  }
}

# Convert results to data frame
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Save results
  output_file <- paste0(mydir, "/fixed_window_results_", mychr, ".RDS")
  saveRDS(results_df, output_file)
  
  cat("\n=== SCANNING COMPLETE ===\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total estimates:", nrow(results_df), "\n")
  cat("Output format: chr | pos | sample | window_size | founder | freq\n")
  
  # Show summary
  cat("\nSummary by window size:\n")
  summary_stats <- results_df %>%
    group_by(window_size) %>%
    summarize(
      n_estimates = n(),
      mean_freq = mean(freq, na.rm = TRUE),
      sd_freq = sd(freq, na.rm = TRUE)
    ) %>%
    arrange(window_size)
  
  print(summary_stats)
  
} else {
  cat("\nNo successful estimates obtained!\n")
}
