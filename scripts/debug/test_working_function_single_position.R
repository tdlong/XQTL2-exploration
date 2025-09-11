#!/usr/bin/env Rscript

# Simple test of the working haplotype estimation function at a single position
# This just tests that the working code works before we try to modify it

suppressPackageStartupMessages({
  library(tidyverse)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/debug/test_working_function_single_position.R <chr> <param_file> <output_dir> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
test_position <- as.numeric(args[4])

cat("=== TESTING WORKING FUNCTION AT SINGLE POSITION ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Test position:", test_position, "\n\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("H cutoff:", h_cutoff, "\n\n")

# Source the working functions
source("scripts/production/haplotype_estimation_functions.R")

# Define euchromatin boundaries (same as the main script)
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

# Load observed data from REFALT files (same as the main script)
cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data file not found: ", refalt_file)
}

# Load REFALT data (same as existing working code)
refalt_data <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(refalt_data), "rows\n")

# Transform to frequencies (same as existing working code)
cat("Converting counts to frequencies...\n")
refalt_processed <- refalt_data %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    total_count = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

# Filter for high-quality SNPs (same as existing working code)
cat("Filtering for high-quality SNPs...\n")
good_snps <- refalt_processed %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(total_count == 0),
    not_fixed = sum(total_count != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Get valid SNPs for euchromatin (same as existing working code)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

# Filter to euchromatin and only the samples we actually processed (same as existing working code)
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Test the working function at the specified position
test_sample <- names_in_bam[1]  # Use first sample
cat("=== TESTING WORKING FUNCTION ===\n")
cat("Test position:", test_position, "\n")
cat("Test sample:", test_sample, "\n\n")

# Call the working function
working_result <- estimate_haplotypes(test_position, test_sample, observed_euchromatic, 
                                    founders, h_cutoff, "adaptive", 
                                    NULL, chr, verbose = 2)

cat("\n=== WORKING FUNCTION RESULTS ===\n")
cat("Estimate OK:", working_result$estimate_OK, "\n")
cat("Final window size:", working_result$final_window_size, "\n")
cat("Number of SNPs:", working_result$n_snps, "\n")
cat("Haplotype frequencies:\n")
if (!is.null(working_result$haplotype_freqs)) {
  print(working_result$haplotype_freqs)
} else {
  cat("NULL\n")
}

cat("\n=== TEST COMPLETE ===\n")



