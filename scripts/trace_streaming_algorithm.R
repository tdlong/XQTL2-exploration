#!/usr/bin/env Rscript

# Trace Streaming Algorithm Step-by-Step
# Follows the exact logic to find the bug

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"
sample_name <- "GJ_3_1"

cat("=== Tracing Streaming Algorithm ===\n")
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

cat("✓ Sample haplotypes loaded:", nrow(sample_haplotypes), "rows\n")
cat("Position range:", min(sample_haplotypes$pos), "-", max(sample_haplotypes$pos), "bp\n\n")

# Load SNP data
refalt_file <- file.path(output_dir, paste0("df3.", chr, ".RDS"))
df2 <- read_rds(refalt_file)

# Filter SNPs to euchromatin
good_snps <- df2 %>%
  filter(name %in% c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end)

cat("✓ Valid euchromatic SNPs:", nrow(valid_snps), "\n")
cat("SNP position range:", min(valid_snps$POS), "-", max(valid_snps$POS), "bp\n\n")

# Convert haplotypes to wide format
haplotype_freqs <- sample_haplotypes %>%
  pivot_wider(names_from = founder, values_from = freq, values_fill = NA)

# Get unique haplotype positions (sorted)
haplotype_positions <- sort(unique(haplotype_freqs$pos))
cat("Haplotype positions:", length(haplotype_positions), "\n")
cat("First 10:", paste(head(haplotype_positions, 10), collapse = ", "), "\n")
cat("Last 10:", paste(tail(haplotype_positions, 10), collapse = ", "), "\n\n")

# Get SNP positions to test
snp_positions <- valid_snps %>% distinct(POS) %>% pull(POS) %>% sort()
cat("SNP positions:", length(snp_positions), "\n")
cat("First 10:", paste(head(snp_positions, 10), collapse = ", "), "\n")
cat("Last 10:", paste(tail(snp_positions, 10), collapse = ", "), "\n\n")

# Test the EXACT streaming algorithm logic
cat("=== Tracing Streaming Algorithm Logic ===\n")

# Initialize results
interpolated_results <- list()

# Initialize interval tracking
current_left_idx <- 1
current_right_idx <- 2

cat("Starting indices: left =", current_left_idx, ", right =", current_right_idx, "\n")
cat("Starting haplotype positions: [", haplotype_positions[current_left_idx], ",", haplotype_positions[current_right_idx], "]\n\n")

# Process first 20 SNPs to trace the logic
for (i in 1:20) {
  snp_pos <- snp_positions[i]
  cat("=== SNP", i, "at position", snp_pos, "===\n")
  
  cat("  Current interval indices: left =", current_left_idx, ", right =", current_right_idx, "\n")
  cat("  Current haplotype positions: [", haplotype_positions[current_left_idx], ",", haplotype_positions[current_right_idx], "]\n")
  
  # Find the haplotype interval containing this SNP
  cat("  Looking for interval containing SNP", snp_pos, "...\n")
  
  while (current_right_idx <= length(haplotype_positions) && 
         haplotype_positions[current_right_idx] < snp_pos) {
    cat("    While loop: haplotype", haplotype_positions[current_right_idx], "< SNP", snp_pos, "\n")
    cat("    Advancing interval...\n")
    current_left_idx <- current_right_idx
    current_right_idx <- current_right_idx + 1
    cat("    New indices: left =", current_left_idx, ", right =", current_right_idx, "\n")
    
    if (current_right_idx <= length(haplotype_positions)) {
      cat("    New haplotype positions: [", haplotype_positions[current_left_idx], ",", haplotype_positions[current_right_idx], "]\n")
    } else {
      cat("    Right index beyond haplotype range!\n")
    }
  }
  
  # Check if we have a valid interval
  if (current_left_idx > length(haplotype_positions) || 
      current_right_idx > length(haplotype_positions)) {
    cat("  ❌ Invalid interval indices: left =", current_left_idx, ", right =", current_right_idx, "\n")
    cat("  Haplotype range: 1 to", length(haplotype_positions), "\n")
    break
  }
  
  left_pos <- haplotype_positions[current_left_idx]
  right_pos <- haplotype_positions[current_right_idx]
  
  cat("  ✓ Final interval: [", left_pos, ",", right_pos, "]\n")
  cat("  SNP", snp_pos, "is between haplotypes\n")
  
  # Check if SNP is actually in the interval
  if (snp_pos >= left_pos && snp_pos <= right_pos) {
    cat("  ✓ SNP position is valid for interpolation\n")
    
    # Get founder columns that exist
    existing_founders <- intersect(c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8"), names(haplotype_freqs))
    cat("  Founders:", paste(existing_founders, collapse = ", "), "\n")
    
    # Extract frequencies for left and right positions
    left_freqs <- haplotype_freqs %>% 
      filter(pos == left_pos) %>% 
      select(all_of(existing_founders))
    
    right_freqs <- haplotype_freqs %>% 
      filter(pos == right_pos) %>% 
      select(all_of(existing_founders))
    
    cat("  Left frequencies:", paste(round(as.numeric(left_freqs[1, ]), 3), collapse = ", "), "\n")
    cat("  Right frequencies:", paste(round(as.numeric(right_freqs[1, ]), 3), collapse = ", "), "\n")
    
    # Convert to numeric vectors
    left_freqs_numeric <- as.numeric(left_freqs[1, ])
    right_freqs_numeric <- as.numeric(right_freqs[1, ])
    
    # If both sides are all NA, skip
    if (all(is.na(left_freqs_numeric)) && all(is.na(right_freqs_numeric))) { 
      cat("  ⚠️  Both sides NA - skipping\n")
      next 
    }
    
    # Interpolation weight
    alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
    cat("  Interpolation weight (alpha):", round(alpha, 3), "\n")
    
    # Vectorized interpolation
    interpolated_freqs <- rep(NA, length(existing_founders))
    
    for (j in seq_along(existing_founders)) {
      left_val <- left_freqs_numeric[j]
      right_val <- right_freqs_numeric[j]
      
      if (is.na(left_val) && is.na(right_val)) {
        interpolated_freqs[j] <- NA
      } else if (is.na(left_val)) {
        interpolated_freqs[j] <- right_val
      } else if (is.na(right_val)) {
        interpolated_freqs[j] <- left_val
      } else {
        interpolated_freqs[j] <- alpha * left_val + (1 - alpha) * right_val
      }
    }
    
    cat("  Interpolated frequencies:", paste(round(interpolated_freqs, 3), collapse = ", "), "\n")
    
    # Store result
    interpolated_results[[as.character(snp_pos)]] <- as.data.frame(t(interpolated_freqs), col.names = existing_founders)
    cat("  ✓ Successfully interpolated and stored SNP", snp_pos, "\n")
    
  } else {
    cat("  ❌ SNP position is NOT in the interval!\n")
    cat("  SNP:", snp_pos, "should be between", left_pos, "and", right_pos, "\n")
  }
  
  cat("\n")
}

cat("=== Trace Complete ===\n")
cat("Total interpolated SNPs:", length(interpolated_results), "\n")
if (length(interpolated_results) > 0) {
  cat("First few results:\n")
  for (i in 1:min(3, length(interpolated_results))) {
    snp_pos <- names(interpolated_results)[i]
    cat("  SNP", snp_pos, ":", paste(round(as.numeric(interpolated_results[[i]][1, ]), 3), collapse = ", "), "\n")
  }
}
