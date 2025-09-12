#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(limSolve)
library(MASS)

# Source the functions without running main execution
source('helpfiles/JUICE_haplotype_parameters.R')

# Define the wide format data loading function
process_refalt_data_wide <- function(refalt_file, founders) {
  cat("Loading RefAlt data in WIDE format from:", refalt_file, "\n")
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Keep data in wide format - no pivot_longer!
  # Extract sample names and create frequency matrices
  sample_names <- unique(str_sub(names(refalt_data)[-c(1,2)], 5))  # Remove CHROM, POS, extract sample names
  
  # Create REF and ALT matrices for each sample
  ref_cols <- names(refalt_data)[str_sub(names(refalt_data), 1, 3) == "REF" & str_sub(names(refalt_data), 5) %in% sample_names]
  alt_cols <- names(refalt_data)[str_sub(names(refalt_data), 1, 3) == "ALT" & str_sub(names(refalt_data), 5) %in% sample_names]
  
  ref_matrix <- as.matrix(refalt_data[, ref_cols])
  alt_matrix <- as.matrix(refalt_data[, alt_cols])
  
  # Calculate frequencies and read counts
  freq_matrix <- ref_matrix / (ref_matrix + alt_matrix)
  n_matrix <- ref_matrix + alt_matrix
  
  # Quality filtering: remove positions with too few reads or extreme frequencies
  valid_positions <- rowSums(n_matrix >= 10) == ncol(n_matrix) &  # All samples have >= 10 reads
                    rowSums(freq_matrix >= 0.03 & freq_matrix <= 0.97, na.rm = TRUE) == ncol(freq_matrix)  # All samples have reasonable frequencies
  
  # Filter to valid positions
  refalt_wide <- refalt_data[valid_positions, ]
  freq_matrix <- freq_matrix[valid_positions, ]
  n_matrix <- n_matrix[valid_positions, ]
  
  # Set column names
  colnames(freq_matrix) <- sample_names
  colnames(n_matrix) <- sample_names
  
  cat("✓ Processed", nrow(refalt_wide), "positions for", length(sample_names), "samples in WIDE format\n")
  
  return(list(
    refalt = refalt_wide,
    freq_matrix = freq_matrix,
    n_matrix = n_matrix,
    sample_names = sample_names
  ))
}

# Test the wide format data loading
cat("=== TESTING WIDE FORMAT DATA LOADING ===\n")
wide_data <- process_refalt_data_wide('process/JUICE/RefAlt.chr2R.txt', founders)

cat('\n=== WIDE DATA STRUCTURE ===\n')
cat('RefAlt dimensions:', dim(wide_data$refalt), '\n')
cat('Freq matrix dimensions:', dim(wide_data$freq_matrix), '\n')
cat('Sample names:', paste(wide_data$sample_names, collapse=', '), '\n')
cat('First few positions:', head(wide_data$refalt$POS, 5), '\n')

# Test a single window
pos <- 5400000
window_start <- max(1, pos - 300000)
window_end <- pos + 300000
cat('\n=== WINDOW TEST ===\n')
cat('Position:', pos, '\n')
cat('Window:', window_start, 'to', window_end, '\n')

positions <- wide_data$refalt$POS
window_positions <- positions >= window_start & positions <= window_end
cat('Positions in window:', sum(window_positions), '\n')

if(sum(window_positions) > 0) {
  window_freq <- wide_data$freq_matrix[window_positions, , drop = FALSE]
  cat('Window freq matrix dimensions:', dim(window_freq), '\n')
  cat('Samples in window data:', colnames(window_freq), '\n')
  cat('Founders present:', sum(founders %in% colnames(window_freq)), 'of', length(founders), '\n')
  cat('Founders:', paste(founders, collapse=', '), '\n')
  cat('Sample names in data:', paste(colnames(window_freq), collapse=', '), '\n')
}
