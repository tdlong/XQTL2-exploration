#!/usr/bin/env Rscript

# Simple script: load df4 data and call est_haps_var directly
# Usage: Rscript scripts/ErrMatrix/simple_estimator_test.R <data_file> <metadata_file> <sample_name>

library(tidyverse)
library(limSolve)
library(MASS)
library(purrr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript scripts/ErrMatrix/simple_estimator_test.R <data_file> <metadata_file> <sample_name>\n")
  cat("Example: Rscript scripts/ErrMatrix/simple_estimator_test.R position_data_chr3R_19610000.RDS position_metadata_chr3R_19610000.RDS Rep01_W_F\n")
  quit(status = 1)
}

data_file <- args[1]
metadata_file <- args[2]
sample_name <- args[3]

cat("=== SIMPLE ESTIMATOR TEST ===\n")
cat("Data file:", data_file, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("Sample:", sample_name, "\n\n")

# Load the data
df4 <- readRDS(data_file)
metadata <- readRDS(metadata_file)

cat("Loaded data:\n")
cat("  SNPs:", nrow(df4), "\n")
cat("  Founders:", paste(metadata$founders, collapse=", "), "\n\n")

# Load BASE_VAR_WIDE.R functions
cat("Loading BASE_VAR_WIDE.R functions...\n")
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Call est_haps_var directly
cat("Calling est_haps_var...\n")
result <- est_haps_var(
  testing_position = metadata$position,
  sample_name = sample_name,
  df3 = df4,
  founders = metadata$founders,
  h_cutoff = 4,
  method = "adaptive",
  chr = metadata$chr,
  verbose = 2
)

cat("\n=== RESULT ===\n")
cat("Groups:", paste(result$Groups, collapse=","), "\n")
cat("Haplotypes:\n")
print(result$Haps)
cat("Error matrix:\n")
print(result$Err)
cat("Names:", paste(result$Names, collapse=","), "\n")

cat("\nDone.\n")
