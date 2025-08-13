#!/usr/bin/env Rscript

# Debug Interpolation Process
# Examine exactly what happens during frequency lookup

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"
sample_name <- "GJ_3_1"

cat("=== Debug Interpolation Process ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n\n")

# Load haplotype results
h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
haplotype_results <- read_rds(haplotype_file) %>%
  filter(h_cutoff == !!h_cutoff)

# Filter to sample and euchromatin
euchromatin_start <- 5398184
euchromatin_end <- 24684540

sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name, 
         pos >= euchromatin_start, 
         pos <= euchromatin_end)

cat("Sample haplotypes in euchromatin:", nrow(sample_haplotypes), "\n")

# Test specific SNP positions from the trace
test_snps <- c(5398203, 5398320, 5398508, 5398594, 5398610, 5398611, 5398614, 5398628, 5398636, 5398645, 5398648, 5398650, 5398818, 5399011, 5399147, 5399818, 5400054, 5400121, 5400151, 5400166)

cat("\n=== Testing SNP Positions ===\n")
for (snp_pos in test_snps) {
  cat("SNP at position", snp_pos, ":\n")
  
  # Find haplotype positions that bracket this SNP
  left_haplo <- sample_haplotypes %>% 
    filter(pos <= snp_pos) %>% 
    arrange(desc(pos)) %>% 
    slice_head(n = 1)
  
  right_haplo <- sample_haplotypes %>% 
    filter(pos >= snp_pos) %>% 
    arrange(pos) %>% 
    slice_head(n = 1)
  
  if (nrow(left_haplo) > 0 && nrow(right_haplo) > 0) {
    cat("  Left haplotype: pos =", left_haplo$pos[1], "\n")
    cat("  Right haplotype: pos =", right_haplo$pos[1], "\n")
    
    # Check if SNP is between haplotypes
    if (left_haplo$pos[1] <= snp_pos && snp_pos <= right_haplo$pos[1]) {
      cat("  ✓ SNP is between haplotypes\n")
      
      # Get frequencies for left haplotype
      left_freqs <- sample_haplotypes %>% 
        filter(pos == left_haplo$pos[1]) %>%
        select(founder, freq) %>%
        arrange(founder)
      
      # Get frequencies for right haplotype  
      right_freqs <- sample_haplotypes %>% 
        filter(pos == right_haplo$pos[1]) %>%
        select(founder, freq) %>%
        arrange(founder)
      
      cat("  Left frequencies:\n")
      for (i in 1:nrow(left_freqs)) {
        cat("    ", left_freqs$founder[i], ":", left_freqs$freq[i], "\n")
      }
      
      cat("  Right frequencies:\n")
      for (i in 1:nrow(right_freqs)) {
        cat("    ", right_freqs$founder[i], ":", right_freqs$freq[i], "\n")
      }
      
      # Check for NA frequencies
      left_na <- sum(is.na(left_freqs$freq))
      right_na <- sum(is.na(right_freqs$freq))
      
      if (left_na > 0 || right_na > 0) {
        cat("  ⚠️  Found NA frequencies: Left =", left_na, "Right =", right_na, "\n")
      } else {
        cat("  ✓ All frequencies are valid\n")
      }
      
    } else {
      cat("  ❌ SNP is NOT between haplotypes\n")
    }
  } else {
    cat("  ❌ Missing haplotype data\n")
  }
  cat("\n")
}

# Check the specific interval that should work
cat("=== Checking Interval [5400000, 5410000] ===\n")
left_pos <- 5400000
right_pos <- 5410000

left_data <- sample_haplotypes %>% filter(pos == left_pos)
right_data <- sample_haplotypes %>% filter(pos == right_pos)

cat("Left position", left_pos, "has", nrow(left_data), "rows\n")
cat("Right position", right_pos, "has", nrow(right_data), "rows\n")

if (nrow(left_data) > 0 && nrow(right_data) > 0) {
  cat("\nLeft frequencies:\n")
  for (i in 1:nrow(left_data)) {
    cat("  ", left_data$founder[i], ":", left_data$freq[i], "\n")
  }
  
  cat("\nRight frequencies:\n")
  for (i in 1:nrow(right_data)) {
    cat("  ", right_data$founder[i], ":", right_data$freq[i], "\n")
  }
  
  # Check for duplicates
  cat("\nChecking for duplicates:\n")
  left_dups <- left_data %>% count(founder) %>% filter(n > 1)
  right_dups <- right_data %>% count(founder) %>% filter(n > 1)
  
  if (nrow(left_dups) > 0) {
    cat("  Left position has duplicate founders:\n")
    print(left_dups)
  }
  
  if (nrow(right_dups) > 0) {
    cat("  Right position has duplicate founders:\n")
    print(right_dups)
  }
}

cat("\n=== Debug Complete ===\n")
