#!/usr/bin/env Rscript

# Debug Interpolation Step by Step
# Examine specific SNPs to see what's going wrong

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "fixed_500kb"
sample_name <- "GJ_3_1"

cat("=== Debug Interpolation Step by Step ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n\n")

# Load data
if (grepl("^fixed_", estimator)) {
  window_size <- as.numeric(gsub("fixed_|kb", "", estimator)) * 1000
  haplotype_file <- file.path(output_dir, paste0("fixed_window_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file) %>%
    filter(window_size == !!window_size)
}

refalt_file <- file.path(output_dir, paste0("df3.", chr, ".RDS"))
df2 <- read_rds(refalt_file)

# Define euchromatin boundaries
euchromatin_start <- 5398184
euchromatin_end <- 24684540

# Get haplotype positions for this sample
sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name)

haplotype_positions <- sort(unique(sample_haplotypes$pos))

# Get observed frequencies for this sample
observed_data <- df2 %>%
  filter(name == sample_name, 
         POS >= euchromatin_start, 
         POS <= euchromatin_end) %>%
  select(POS, freq) %>%
  rename(observed = freq)

# Sample a few SNPs for detailed debugging
set.seed(123)
test_snps <- sample(observed_data$POS, 5)

cat("=== Testing 5 Random SNPs ===\n")
for (i in seq_along(test_snps)) {
  snp_pos <- test_snps[i]
  cat("\n--- SNP", i, "at position", snp_pos, "---\n")
  
  # Get observed frequency
  observed_freq <- observed_data %>% filter(POS == snp_pos) %>% pull(observed)
  cat("Observed frequency:", observed_freq, "\n")
  
  # Find flanking haplotypes
  left_haplo <- sample_haplotypes %>% 
    filter(pos <= snp_pos) %>% 
    arrange(desc(pos)) %>% 
    slice_head(n = 1)
  
  right_haplo <- sample_haplotypes %>% 
    filter(pos >= snp_pos) %>% 
    arrange(pos) %>% 
    slice_head(n = 1)
  
  if (nrow(left_haplo) > 0 && nrow(right_haplo) > 0) {
    cat("Left haplotype position:", left_haplo$pos[1], "\n")
    cat("Right haplotype position:", right_haplo$pos[1], "\n")
    
    # Get founder frequencies for left haplotype
    left_freqs <- sample_haplotypes %>% 
      filter(pos == left_haplo$pos[1]) %>%
      select(founder, freq) %>%
      arrange(founder)
    
    cat("Left haplotype founder frequencies:\n")
    for (j in 1:nrow(left_freqs)) {
      cat("  ", left_freqs$founder[j], ":", left_freqs$freq[j], "\n")
    }
    
    # Get founder frequencies for right haplotype
    right_freqs <- sample_haplotypes %>% 
      filter(pos == right_haplo$pos[1]) %>%
      select(founder, freq) %>%
      arrange(founder)
    
    cat("Right haplotype founder frequencies:\n")
    for (j in 1:nrow(right_freqs)) {
      cat("  ", right_freqs$founder[j], ":", right_freqs$freq[j], "\n")
    }
    
    # Calculate interpolation
    alpha <- (right_haplo$pos[1] - snp_pos) / (right_haplo$pos[1] - left_haplo$pos[1])
    cat("Interpolation weight (alpha):", round(alpha, 3), "\n")
    
    # Interpolate each founder
    cat("Interpolated founder frequencies:\n")
    for (founder in unique(c(left_freqs$founder, right_freqs$founder))) {
      left_val <- left_freqs %>% filter(founder == !!founder) %>% pull(freq)
      right_val <- right_freqs %>% filter(founder == !!founder) %>% pull(freq)
      
      if (length(left_val) == 0) left_val <- NA
      if (length(right_val) == 0) right_val <- NA
      
      if (is.na(left_val) && is.na(right_val)) {
        interpolated <- NA
      } else if (is.na(left_val)) {
        interpolated <- right_val
      } else if (is.na(right_val)) {
        interpolated <- left_val
      } else {
        interpolated <- alpha * left_val + (1 - alpha) * right_val
      }
      
      cat("  ", founder, ":", round(interpolated, 4), 
          " (left:", round(left_val, 4), "right:", round(right_val, 4), ")\n")
    }
    
    # Calculate imputed frequency (simple average for now)
    interpolated_freqs <- left_freqs %>%
      left_join(right_freqs, by = "founder", suffix = c("_left", "_right")) %>%
      mutate(
        interpolated = case_when(
          is.na(freq_left) & is.na(freq_right) ~ NA_real_,
          is.na(freq_left) ~ freq_right,
          is.na(freq_right) ~ freq_left,
          TRUE ~ alpha * freq_left + (1 - alpha) * freq_right
        )
      )
    
    # For now, just use first founder as example
    if (nrow(interpolated_freqs) > 0) {
      first_founder_interp <- interpolated_freqs$interpolated[1]
      cat("First founder interpolated frequency:", round(first_founder_interp, 4), "\n")
      cat("Error vs observed:", round(abs(first_founder_interp - observed_freq), 4), "\n")
    }
    
  } else {
    cat("❌ No flanking haplotypes found\n")
  }
}

cat("\n=== Summary ===\n")
cat("This will show us exactly what frequencies are being used\n")
cat("and where the interpolation is going wrong.\n")
