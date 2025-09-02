#!/usr/bin/env Rscript

# SNP Imputation Coverage Diagnostic
# This script analyzes haplotype results and REFALT data to understand
# why some SNPs get imputed and others get NA values
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
cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load haplotype results
cat("Loading haplotype results...\n")
if (grepl("^fixed_", estimator)) {
  window_size_kb <- as.numeric(gsub("fixed_", "", gsub("kb", "", estimator)))
  haplotype_file <- file.path(output_dir, paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
} else if (grepl("^adaptive_h", estimator)) {
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
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
refalt_file <- file.path(dirname(output_dir), paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT file not found:", refalt_file)
}

# Read REFALT data in chunks for memory efficiency
refalt_cols <- c("CHROM", "POS", "REF", "ALT", names_in_bam)
refalt_data <- read_tsv(refalt_file, col_names = refalt_cols, col_types = cols(), show_col_types = FALSE)
cat("✓ REFALT data loaded:", nrow(refalt_data), "rows\n")

# Filter to euchromatin and high-quality SNPs
cat("Filtering to euchromatin SNPs...\n")
euchromatic_snps <- refalt_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(!is.na(REF), !is.na(ALT), REF != "", ALT != "") %>%
  filter(REF != ALT)  # Remove monomorphic SNPs

cat("✓ Euchromatic SNPs:", nrow(euchromatic_snps), "\n\n")

# Analyze haplotype coverage for each sample
cat("=== Haplotype Coverage Analysis ===\n")

for (sample_name in names_in_bam) {
  cat("\n--- Sample:", sample_name, "---\n")
  
  # Get haplotype results for this sample
  sample_haplotypes <- haplotype_results %>%
    filter(sample == sample_name) %>%
    filter(pos >= euchromatin_start, pos <= euchromatin_end)
  
  cat("Haplotype positions:", nrow(sample_haplotypes), "\n")
  
  # Check for NA haplotypes
  na_positions <- sample_haplotypes %>%
    filter(if_any(all_of(founders), is.na))
  
  cat("Positions with NA haplotypes:", nrow(na_positions), "\n")
  
  if (nrow(na_positions) > 0) {
    # Get unique positions with NA haplotypes
    na_pos_list <- unique(na_positions$pos)
    cat("First 5 NA positions:", paste(head(na_pos_list, 5), collapse = ", "), "\n")
    
    # Count SNPs in intervals flanked by NA haplotypes
    snps_in_na_intervals <- 0
    na_intervals <- 0
    
    # Sort haplotype positions
    haplotype_positions <- sort(sample_haplotypes$pos)
    
    # Find intervals between haplotype positions
    for (i in 1:(length(haplotype_positions) - 1)) {
      left_pos <- haplotype_positions[i]
      right_pos <- haplotype_positions[i + 1]
      
      # Check if both positions have NA haplotypes
      left_has_na <- any(is.na(sample_haplotypes %>% filter(pos == left_pos) %>% select(all_of(founders))))
      right_has_na <- any(is.na(sample_haplotypes %>% filter(pos == right_pos) %>% select(all_of(founders))))
      
      if (left_has_na && right_has_na) {
        na_intervals <- na_intervals + 1
        
        # Count SNPs in this interval
        snps_in_interval <- euchromatic_snps %>%
          filter(POS > left_pos, POS < right_pos) %>%
          nrow()
        
        snps_in_na_intervals <- snps_in_na_intervals + snps_in_interval
        
        if (na_intervals <= 5) {
          cat("  NA interval", na_intervals, ":", left_pos, "-", right_pos, "(", snps_in_interval, "SNPs)\n")
        }
      }
    }
    
    cat("Total NA intervals:", na_intervals, "\n")
    cat("SNPs in NA intervals:", snps_in_na_intervals, "\n")
  }
  
  # Check founder state coverage
  founder_states <- euchromatic_snps %>%
    filter(name %in% founders) %>%
    select(CHROM, POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq) %>%
    arrange(POS)
  
  missing_founder_states <- founder_states %>%
    filter(if_any(all_of(founders), is.na))
  
  cat("SNPs with missing founder states:", nrow(missing_founder_states), "\n")
  
  # Summary for this sample
  total_snps <- nrow(euchromatic_snps)
  snps_with_haplotypes <- nrow(sample_haplotypes)
  snps_with_founder_states <- nrow(founder_states) - nrow(missing_founder_states)
  
  cat("Coverage summary:\n")
  cat("  Total euchromatic SNPs:", total_snps, "\n")
  cat("  SNPs with haplotypes:", snps_with_haplotypes, "\n")
  cat("  SNPs with founder states:", snps_with_founder_states, "\n")
  cat("  SNPs in NA intervals:", snps_in_na_intervals, "\n")
  cat("  Expected imputations:", min(snps_with_haplotypes, snps_with_founder_states) - snps_in_na_intervals, "\n")
}

cat("\n=== Diagnostic Complete ===\n")
cat("This analysis shows why some SNPs get imputed and others don't.\n")
cat("The key is understanding haplotype coverage and NA intervals.\n")
