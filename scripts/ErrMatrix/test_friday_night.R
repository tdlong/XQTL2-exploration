#!/usr/bin/env Rscript

# Test script to reproduce Friday night results
# This runs BASE_VAR_WIDE.R directly like the Friday night run did

# Set up the same arguments as Friday night
chr <- "chr3R"
method <- "adaptive"
parameter <- 4
output_dir <- "process/ZINC2"
param_file <- "helpfiles/ZINC2_haplotype_parameters.R"

cat("=== TESTING FRIDAY NIGHT REPRODUCTION ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Run BASE_VAR_WIDE.R with the same arguments as Friday night
# but in debug mode to limit to one position
cat("Running BASE_VAR_WIDE.R with Friday night arguments...\n")

# Set up command line arguments for BASE_VAR_WIDE.R
args <- c(chr, method, parameter, output_dir, param_file, "--debug", "--nonverbose")

# Override commandArgs to simulate the Friday night call
commandArgs <- function(trailingOnly = TRUE) args

# Source and run BASE_VAR_WIDE.R
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")

cat("\nDone.\n")
