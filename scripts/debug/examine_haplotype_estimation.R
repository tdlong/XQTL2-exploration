#!/usr/bin/env Rscript

# Examine Specific Haplotype Estimation
# This script examines a specific haplotype position to see what SNPs are being used
#
# Usage: Rscript examine_haplotype_estimation.R <chr> <param_file> <output_dir> <estimator> <position>

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript examine_haplotype_estimation.R <chr> <param_file> <output_dir> <estimator> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]
target_position <- as.numeric(args[5])

cat("=== Examine Specific Haplotype Estimation ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Estimator:", estimator, "\n")
cat("Target position:", target_position, "\n\n")

# Load parameter file
cat("Loading parameter file...\n")
source(param_file)
cat("✓ Parameter file loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("Samples to process:", paste(names_in_bam, collapse = ", "), "\n\n")

# Define euchromatin boundaries
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]

# Load haplotype results
cat("Loading haplotype results...\n")
if (grepl("^fixed_", estimator)) {
  window_size_kb <- as.numeric(gsub("fixed_", "", gsub("kb", "", estimator)))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
  cat("Fixed window size:", window_size_kb, "kb\n")
} else if (grepl("^adaptive_h", estimator)) {
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
  cat("Adaptive h_cutoff:", h_cutoff, "\n")
} else {
  stop("Invalid estimator format")
}

if (!file.exists(haplotype_file)) {
  stop("Haplotype file not found:", haplotype_file)
}

haplotype_results <- read_rds(haplotype_file)
cat("✓ Haplotype results loaded:", nrow(haplotype_results), "rows\n")

# Load REFALT data
cat("Loading REFALT data...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT file not found:", refalt_file)
}

# Load and process REFALT data
df <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(df), "rows\n")

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

cat("✓ Processed REFALT data:", nrow(df2), "rows\n")

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

# Get valid SNPs for evaluation (euchromatin only)
valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(df2, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n\n")

# Examine the specific haplotype position
cat("=== Examining Haplotype at Position", target_position, "===\n")

# Get haplotype results for this position
haplotype_at_position <- haplotype_results %>%
  filter(pos == target_position)

if (nrow(haplotype_at_position) == 0) {
  cat("❌ No haplotype results found at position", target_position, "\n")
} else {
  cat("✓ Haplotype results found at position", target_position, "\n")
  
  # Show haplotype frequencies for each sample
  for (sample_name in names_in_bam) {
    cat("\n--- Sample:", sample_name, "---\n")
    
    sample_haplotype <- haplotype_at_position %>%
      filter(sample == sample_name)
    
    if (nrow(sample_haplotype) == 0) {
      cat("  No haplotype data for this sample\n")
    } else {
      # Check for NA haplotypes
      has_na <- any(is.na(sample_haplotype %>% select(all_of(founders))))
      cat("  Has NA haplotypes:", has_na, "\n")
      
      if (has_na) {
        cat("  NA founder columns:", paste(names(which(is.na(sample_haplotype %>% select(all_of(founders)) %>% as.numeric()))), collapse = ", "), "\n")
      }
      
      # Show haplotype frequencies
      cat("  Haplotype frequencies:\n")
      for (founder in founders) {
        if (founder %in% names(sample_haplotype)) {
          freq <- sample_haplotype %>% pull(founder)
          cat("    ", founder, ":", freq, "\n")
        } else {
          cat("    ", founder, ": column missing\n")
        }
      }
    }
  }
}

# Now examine the SNP coverage around this position
cat("\n=== SNP Coverage Analysis ===\n")

# Calculate expected window boundaries for fixed window
if (grepl("^fixed_", estimator)) {
  window_size_bp <- window_size_kb * 1000
  expected_left <- target_position - (window_size_bp / 2)
  expected_right <- target_position + (window_size_bp / 2)
  
  cat("Expected 20kb window boundaries:\n")
  cat("  Left boundary:", expected_left, "\n")
  cat("  Right boundary:", expected_right, "\n")
  cat("  Window size:", window_size_bp, "bp\n\n")
  
  # Count SNPs in the expected window
  snps_in_window <- valid_snps %>%
    filter(POS >= expected_left, POS <= expected_right) %>%
    distinct(CHROM, POS) %>%
    nrow()
  
  cat("SNPs in expected", window_size_kb, "kb window:", snps_in_window, "\n")
  
  # Show SNP distribution
  snp_distribution <- valid_snps %>%
    filter(POS >= expected_left, POS <= expected_right) %>%
    distinct(CHROM, POS) %>%
    mutate(
      distance_from_center = abs(POS - target_position)
    ) %>%
    arrange(distance_from_center)
  
  cat("SNP distribution around position", target_position, ":\n")
  cat("  Closest SNP:", min(snp_distribution$distance_from_center), "bp away\n")
  cat("  Farthest SNP:", max(snp_distribution$distance_from_center), "bp away\n")
  cat("  SNPs within 10kb:", sum(snp_distribution$distance_from_center <= 10000), "\n")
  cat("  SNPs within 20kb:", sum(snp_distribution$distance_from_center <= 20000), "\n")
  
  # Show first few SNPs
  cat("\nFirst 10 SNPs in window:\n")
  print(head(snp_distribution, 10))
  
} else {
  cat("Adaptive window - cannot calculate fixed boundaries\n")
}

cat("\n=== Analysis Complete ===\n")
cat("This shows what SNPs the haplotype estimator actually used.\n")
cat("Check if the window boundaries match the expected", window_size_kb, "kb size.\n")
