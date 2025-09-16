#!/usr/bin/env Rscript

# Extract single position data using EXACT logic from BASE_VAR_WIDE.R
# Usage: Rscript scripts/ErrMatrix/extract_single_position_data.R <chr> <position> <sample_name> <output_dir> <param_file>

library(tidyverse)
library(limSolve)
library(MASS)
library(purrr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript scripts/ErrMatrix/extract_single_position_data.R <chr> <position> <sample_name> <output_dir> <param_file>\n")
  cat("Example: Rscript scripts/ErrMatrix/extract_single_position_data.R chr3R 19610000 Rep01_W_F process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R\n")
  quit(status = 1)
}

chr <- args[1]
testing_position <- as.numeric(args[2])
sample_name <- args[3]
output_dir <- args[4]
param_file <- args[5]

cat("=== EXTRACTING SINGLE POSITION DATA ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", testing_position, "\n")
cat("Sample:", sample_name, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Source the Friday night working version of BASE_VAR_WIDE.R
old_interactive <- interactive
interactive <- function() TRUE
# Use the current working version with dplyr fix
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Load parameters (EXACT from BASE_VAR_WIDE.R)
source(param_file)

# Load RefAlt data (EXACT from BASE_VAR_WIDE.R)
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

df3 <- process_refalt_data(refalt_file, founders)

# Pre-subset df3 to df4 for this position (EXACT from BASE_VAR_WIDE.R)
max_window <- 500000  # Largest window size in est_haps_var
window_start <- max(1, testing_position - max_window/2)
window_end <- testing_position + max_window/2

df4 <- df3 %>%
  filter(POS >= window_start & POS <= window_end)

cat("Subsetted to", nrow(df4), "SNPs in window", window_start, "-", window_end, "\n")

# Save the df4 data and parameters (EXACT what gets passed to estimate_haplotypes_list_format)
dump_payload <- list(
  testing_position = testing_position,
  sample_name = sample_name,
  df4 = df4,  # Pre-subsetted data
  founders = founders,
  h_cutoff = h_cutoff,
  method = "adaptive",
  chr = chr
)

# Save to file
dump_file <- paste0("position_data_", chr, "_", testing_position, "_", sample_name, ".RDS")
saveRDS(dump_payload, dump_file)

cat("Saved data to:", dump_file, "\n")
cat("  Position:", testing_position, "\n")
cat("  Sample:", sample_name, "\n")
cat("  SNPs in window:", nrow(df4), "\n")
cat("  Founders:", paste(founders, collapse=", "), "\n")
cat("  h_cutoff:", h_cutoff, "\n")

cat("\nDone.\n")
