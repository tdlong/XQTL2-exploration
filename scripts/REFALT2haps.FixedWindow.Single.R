#!/usr/bin/env Rscript

# REFALT2haps Fixed Window - Single Parameter Version
# This script runs haplotype estimation for a single fixed window size
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
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Fixed Window - Single Parameter ===\n")
cat("Chromosome:", mychr, "\n")
cat("Window size:", window_size_kb, "kb (", window_size_bp, "bp)\n")
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

# Get all non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("Testing window size:", window_size_kb, "kb\n")
cat("Using fixed h_cutoff:", h_cutoff, "\n")
cat("Non-founder samples:", length(non_founder_samples), "\n")
cat("Samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Define scanning positions (500kb to end-500kb, 10kb steps)
chromosome_length <- max(df$POS)
scan_start <- 500000
scan_end <- chromosome_length - 500000
scan_positions <- seq(scan_start, scan_end, by = 10000)

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
  window_start <- test_pos - window_size_bp
  window_end <- test_pos + window_size_bp
  
  # Get SNPs in window
  window_snps <- df3 %>%
    filter(POS >= window_start & POS <= window_end)
  
  if (nrow(window_snps) < 10) {
    # Skip windows with too few SNPs
    next
  }
  
  # Process each sample
  for (sample_name in non_founder_samples) {
    # Get sample data
    sample_data <- window_snps %>%
      filter(name == sample_name) %>%
      select(POS, freq, N)
    
    if (nrow(sample_data) < 5) next  # Skip if too few SNPs
    
    # Get founder data for this window
    founder_data <- window_snps %>%
      filter(name %in% founders) %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (ncol(founder_data) < length(founders) + 2) next  # Skip if missing founders
    
    # Prepare founder matrix (exclude POS column)
    founder_matrix <- founder_data %>%
      select(-POS) %>%
      as.matrix()
    
    # Remove rows with any NA values
    complete_rows <- complete.cases(founder_matrix)
    if (sum(complete_rows) < 5) next  # Skip if too few complete rows
    
    founder_matrix <- founder_matrix[complete_rows, ]
    sample_freqs <- sample_data$freq[match(founder_data$POS[complete_rows], sample_data$POS)]
    
    # Remove NA sample frequencies
    valid_indices <- !is.na(sample_freqs)
    if (sum(valid_indices) < 5) next  # Skip if too few valid frequencies
    
    founder_matrix <- founder_matrix[valid_indices, ]
    sample_freqs <- sample_freqs[valid_indices]
    
    # Solve constrained least squares
    tryCatch({
      result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, 
                              E = matrix(1, nrow = 1, ncol = ncol(founder_matrix)), 
                              F = 1, G = diag(ncol(founder_matrix)), H = rep(0, ncol(founder_matrix)))
      
      if (result$IsError == 0) {
        # Store results
        result_row <- list(
          chr = mychr,
          pos = test_pos,
          sample = sample_name,
          window_size = window_size_bp,
          n_snps = nrow(window_snps),
          founder_frequencies = result$X
        )
        
        # Add founder frequencies as named columns
        for (i in seq_along(founders)) {
          result_row[[founders[i]]] <- result$X[i]
        }
        
        results_list[[length(results_list) + 1]] <- result_row
      }
    }, error = function(e) {
      # Skip this window if there's an error
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
  cat("Total haplotype estimates:", nrow(results_df), "\n")
  cat("Output file:", output_file, "\n")
  cat("File size:", file.size(output_file), "bytes\n")
  
  # Show sample summary
  cat("\nEstimates per sample:\n")
  sample_counts <- results_df %>%
    group_by(sample) %>%
    summarize(n_estimates = n()) %>%
    arrange(desc(n_estimates))
  print(sample_counts)
  
} else {
  cat("\n❌ No haplotype estimates obtained!\n")
  cat("This could be due to:\n")
  cat("- Insufficient SNP coverage in windows\n")
  cat("- Window size too small for reliable estimation\n")
  cat("- Data quality issues\n")
}
