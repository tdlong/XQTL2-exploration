#!/usr/bin/env Rscript

# SNP Imputation Coverage Diagnostic - SIMPLIFIED VERSION
# Reuses existing pipeline code instead of reinventing the wheel
#
# Usage: Rscript diagnose_snp_imputation_coverage.R <chr> <param_file> <output_dir> <estimator>

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript diagnose_snp_imputation_coverage.R <chr> <param_file> <output_dir> <estimator>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]

cat("=== SNP Imputation Coverage Diagnostic ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Estimator:", estimator, "\n\n")

# Load parameter file
cat("Loading parameter file...\n")
source(param_file)
cat("✓ Parameter file loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("Samples to process:", paste(names_in_bam, collapse = ", "), "\n\n")

# Define euchromatin boundaries (copied from existing pipeline)
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load haplotype results (copied from existing pipeline)
cat("Loading haplotype results...\n")
if (grepl("^fixed_", estimator)) {
  window_size_kb <- as.numeric(gsub("fixed_", "", gsub("kb", "", estimator)))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
} else if (grepl("^adaptive_h", estimator)) {
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
} else {
  stop("Invalid estimator format")
}

if (!file.exists(haplotype_file)) {
  stop("Haplotype file not found:", haplotype_file)
}

haplotype_results <- read_rds(haplotype_file)
cat("✓ Haplotype results loaded:", nrow(haplotype_results), "rows\n")

# Load REFALT data (copied from existing pipeline)
cat("Loading REFALT data...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT file not found:", refalt_file)
}

# Load and process REFALT data (copied from existing pipeline)
df <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(df), "rows\n")

# Transform REF/ALT counts to frequencies (copied from existing pipeline)
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

# Filter for high-quality SNPs (copied from existing pipeline)
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

# Now analyze what the existing pipeline would do
cat("=== Analysis of Existing Pipeline Logic ===\n")

for (sample_name in names_in_bam) {
  cat("\n--- Sample:", sample_name, "---\n")
  
  # Get haplotype results for this sample (copied from existing pipeline)
  sample_haplotypes <- haplotype_results %>%
    filter(sample == sample_name) %>%
    filter(pos >= euchromatin_start, pos <= euchromatin_end)
  
  cat("Haplotype positions (every 10kb):", nrow(sample_haplotypes), "\n")
  
  # Check for NA haplotypes
  na_positions <- sample_haplotypes %>%
    filter(if_any(all_of(founders), is.na))
  
  cat("Haplotype positions with NA values:", nrow(na_positions), "\n")
  
  # Count SNPs that can be imputed vs. skipped
  total_snps <- nrow(valid_snps %>% distinct(CHROM, POS))
  snps_imputable <- 0
  snps_skipped_due_to_na <- 0
  na_intervals <- 0
  
  if (nrow(sample_haplotypes) > 1) {
    # Sort haplotype positions
    haplotype_positions <- sort(sample_haplotypes$pos)
    
    # Analyze each interval between haplotype positions
    for (i in 1:(length(haplotype_positions) - 1)) {
      left_pos <- haplotype_positions[i]
      right_pos <- haplotype_positions[i + 1]
      
      # Check if both positions have NA haplotypes
      left_has_na <- any(is.na(sample_haplotypes %>% filter(pos == left_pos) %>% select(all_of(founders))))
      right_has_na <- any(is.na(sample_haplotypes %>% filter(pos == right_pos) %>% select(all_of(founders))))
      
      # Count SNPs in this interval
      snps_in_interval <- valid_snps %>%
        filter(POS > left_pos, POS < right_pos) %>%
        distinct(CHROM, POS) %>%
        nrow()
      
      if (left_has_na && right_has_na) {
        # Both flanking haplotypes are NA - SNPs in this interval get skipped
        na_intervals <- na_intervals + 1
        snps_skipped_due_to_na <- snps_skipped_due_to_na + snps_in_interval
        
        if (na_intervals <= 3) {
          cat("  NA interval", na_intervals, ":", left_pos, "-", right_pos, "(", snps_in_interval, "SNPs skipped)\n")
        }
      } else {
        # At least one flanking haplotype is valid - SNPs in this interval can be imputed
        snps_imputable <- snps_imputable + snps_in_interval
      }
    }
  }
  
  # Check founder state coverage (copied from existing pipeline logic)
  founder_states_wide <- valid_snps %>%
    filter(name %in% founders) %>%
    select(CHROM, POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq) %>%
    arrange(POS)
  
  missing_founder_states <- founder_states_wide %>%
    filter(if_any(all_of(founders), is.na))
  
  cat("SNPs with missing founder states:", nrow(missing_founder_states), "\n")
  
  # Clear summary
  cat("Coverage summary:\n")
  cat("  Total euchromatic SNPs:", total_snps, "\n")
  cat("  Haplotype positions (every 10kb):", nrow(sample_haplotypes), "\n")
  cat("  SNPs that can be imputed:", snps_imputable, "\n")
  cat("  SNPs that get skipped due to NA haplotypes:", snps_skipped_due_to_na, "\n")
  cat("  NA intervals (both sides NA):", na_intervals, "\n")
  cat("  SNPs with founder states:", nrow(founder_states_wide) - nrow(missing_founder_states), "\n")
  cat("  Expected successful imputations:", min(snps_imputable, nrow(founder_states_wide) - nrow(missing_founder_states)), "\n")
}

cat("\n=== Diagnostic Complete ===\n")
cat("This shows what the existing pipeline logic would do.\n")
cat("Key insight: SNPs flanked by NA haplotypes get skipped entirely.\n")
cat("Smaller fixed windows have more NA haplotypes, so more SNPs get skipped.\n")
