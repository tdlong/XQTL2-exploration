#!/usr/bin/env Rscript

# Load extracted position data and run the estimator with verbose output
# Usage: Rscript scripts/ErrMatrix/run_estimator_on_position.R <chr> <position> <sample> <data_file> <metadata_file>

library(tidyverse)
library(limSolve)
library(MASS)
library(purrr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript scripts/ErrMatrix/run_estimator_on_position.R <chr> <position> <sample> <data_file> <metadata_file>\n")
  cat("Example: Rscript scripts/ErrMatrix/run_estimator_on_position.R chr3R 19610000 Rep01_W_F position_data_chr3R_19610000.RDS position_metadata_chr3R_19610000.RDS\n")
  quit(status = 1)
}

chr <- args[1]
position <- as.numeric(args[2])
sample_name <- args[3]
data_file <- args[4]
metadata_file <- args[5]

cat("=== RUNNING ESTIMATOR ON EXTRACTED POSITION DATA ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", position, "\n")
cat("Sample:", sample_name, "\n")
cat("Data file:", data_file, "\n")
cat("Metadata file:", metadata_file, "\n\n")

# Load the extracted data
if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

cat("Loading extracted data...\n")
df4 <- readRDS(data_file)
metadata <- readRDS(metadata_file)

cat("Loaded data:\n")
cat("  SNPs in window:", nrow(df4), "\n")
cat("  Position range:", metadata$pos_range[1], "to", metadata$pos_range[2], "\n")
cat("  Founders:", paste(metadata$founders, collapse=", "), "\n")
cat("  Samples:", paste(metadata$samples, collapse=", "), "\n\n")

# Check if sample exists in data
if (!sample_name %in% metadata$samples) {
  stop("Sample ", sample_name, " not found in data. Available samples: ", paste(metadata$samples, collapse=", "))
}

# Source the Friday night working code
cat("Loading estimator functions...\n")
source("scripts/debug/haplotype_estimation_functions_working_copy.R")

# Set parameters (same as Friday night)
parameter <- 4
method <- "adaptive"

cat("Running estimator with parameters:\n")
cat("  Parameter (h_cutoff):", parameter, "\n")
cat("  Method:", method, "\n")
cat("  Verbose: 2 (detailed output)\n\n")

# Run the estimator
cat("=== ESTIMATOR OUTPUT ===\n")
result <- estimate_haplotypes_list_format(
  pos = position,
  sample_name = sample_name,
  df3 = df4,
  founders = metadata$founders,
  h_cutoff = parameter,
  method = method,
  chr = chr,
  verbose = 2
)

cat("\n=== FINAL RESULT ===\n")
cat("Groups:", paste(result$Groups, collapse=","), "\n")
cat("Haplotypes:\n")
print(result$Haps)
cat("Error matrix:\n")
print(result$Err)
cat("Names:", paste(result$Names, collapse=","), "\n")

# Save the result
result_file <- paste0("estimator_result_", chr, "_", position, "_", sample_name, ".RDS")
saveRDS(result, result_file)
cat("\nSaved result to:", result_file, "\n")

cat("\nDone.\n")
