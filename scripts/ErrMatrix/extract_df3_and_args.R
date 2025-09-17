#!/usr/bin/env Rscript

# Extract df3 and estimator arguments exactly as used by est_haps_var()
# - Builds df3 via process_refalt_data() from BASE_VAR_WIDE.R to ensure parity
# - Saves a single RDS containing: df3 and a list of function args
#
# Usage (on cluster):
# Rscript scripts/ErrMatrix/extract_df3_and_args.R <chr> <param_file> <output_dir> <position> <sample_name> <method> <window_size_bp>
# Example:
# Rscript scripts/ErrMatrix/extract_df3_and_args.R chr3R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 19780000 Rep01_W_F adaptive 100000

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript scripts/ErrMatrix/extract_df3_and_args.R <chr> <param_file> <output_dir> <position> <sample_name> <method> <window_size_bp>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
testing_position <- as.numeric(args[4])
sample_name <- args[5]
method <- args[6]
window_size_bp <- as.numeric(args[7])

cat("=== Extract df3 and arguments for est_haps_var() ===\n")
cat("chr:", chr, "\n")
cat("param_file:", param_file, "\n")
cat("output_dir:", output_dir, "\n")
cat("position:", testing_position, "\n")
cat("sample:", sample_name, "\n")
cat("method:", method, "\n")
cat("window_size_bp:", window_size_bp, "\n\n")

# Load analysis parameters (founders, h_cutoff, etc.)
source(param_file)
if (!exists("founders")) stop("Param file must define 'founders'")
if (!exists("h_cutoff")) {
  h_cutoff <- 4  # Default value used throughout codebase
  cat("h_cutoff not found in param file, using default:", h_cutoff, "\n")
}

# Use the exact df3 constructor from BASE_VAR_WIDE.R
base_script <- "scripts/ErrMatrix/BASE_VAR_WIDE.R"
if (!file.exists(base_script)) stop("Cannot find ", base_script)
source(base_script)

refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) stop("RefAlt file not found: ", refalt_file)

cat("Building df3 via process_refalt_data()...\n")
df3 <- process_refalt_data(refalt_file, founders)
stopifnot("POS" %in% names(df3))
cat("✓ df3 ready:", nrow(df3), "rows x", ncol(df3), "cols\n")

# Package everything needed to call est_haps_var() later
payload <- list(
  df3 = df3,
  args = list(
    testing_position = testing_position,
    sample_name = sample_name,
    founders = founders,
    h_cutoff = h_cutoff,
    method = method,
    window_size_bp = window_size_bp,
    chr = chr,
    verbose = 2
  )
)

outfile <- sprintf("estimator_df3_and_args_%s_%d_%s_%s_win%d.rds", chr, testing_position, sample_name, method, window_size_bp)
saveRDS(payload, outfile)
cat("Saved:", outfile, "\n")


