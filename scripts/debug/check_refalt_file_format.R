#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

# =============================================================================
# Check REFALT File Format
# =============================================================================
# 
# This script checks the format of the REFALT file to understand why
# the parsing is failing.
#
# USAGE:
# Rscript scripts/debug/check_refalt_file_format.R <chr> <output_dir>
#
# EXAMPLE:
# Rscript scripts/debug/check_refalt_file_format.R chr2R process/ZINC2
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/debug/check_refalt_file_format.R <chr> <output_dir>")
}

chr <- args[1]
output_dir <- args[2]

cat("=== CHECKING REFALT FILE FORMAT ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Check if the file exists
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found:", refalt_file)
}

cat("✓ REFALT file found:", refalt_file, "\n")

# Check file size
file_size <- file.info(refalt_file)$size
cat("File size:", round(file_size / 1024^2, 2), "MB\n\n")

# Look at first few lines
cat("=== FIRST 10 LINES ===\n")
first_lines <- readLines(refalt_file, n = 10)
for (i in 1:length(first_lines)) {
  cat(sprintf("%2d: %s\n", i, first_lines[i]))
}

cat("\n=== TRYING DIFFERENT PARSING APPROACHES ===\n")

# Try reading with different column specifications
cat("1. Trying with default column types...\n")
tryCatch({
  refalt_default <- read_tsv(refalt_file, col_names = c("chr", "pos", "name", "ref_count", "alt_count"), 
                            show_col_types = FALSE)
  cat("✓ Default parsing successful:", nrow(refalt_default), "rows\n")
  cat("Column types:\n")
  str(refalt_default)
}, error = function(e) {
  cat("❌ Default parsing failed:", e$message, "\n")
})

cat("\n2. Trying with explicit column types...\n")
tryCatch({
  refalt_explicit <- read_tsv(refalt_file, col_names = c("chr", "pos", "name", "ref_count", "alt_count"), 
                             col_types = "cdcdd", show_col_types = FALSE)
  cat("✓ Explicit parsing successful:", nrow(refalt_explicit), "rows\n")
  cat("Column types:\n")
  str(refalt_explicit)
}, error = function(e) {
  cat("❌ Explicit parsing failed:", e$message, "\n")
})

cat("\n3. Trying with character types first...\n")
tryCatch({
  refalt_char <- read_tsv(refalt_file, col_names = c("chr", "pos", "name", "ref_count", "alt_count"), 
                         col_types = "ccccc", show_col_types = FALSE)
  cat("✓ Character parsing successful:", nrow(refalt_char), "rows\n")
  cat("First few rows:\n")
  print(head(refalt_char, 3))
  
  # Check for parsing issues
  parsing_issues <- problems(refalt_char)
  if (nrow(parsing_issues) > 0) {
    cat("⚠️  Parsing issues found:", nrow(parsing_issues), "rows\n")
    cat("First few issues:\n")
    print(head(parsing_issues, 5))
  }
}, error = function(e) {
  cat("❌ Character parsing failed:", e$message, "\n")
})

cat("\n=== CHECKING SAMPLE NAMES ===\n")
# Try to get unique names
tryCatch({
  refalt_char <- read_tsv(refalt_file, col_names = c("chr", "pos", "name", "ref_count", "alt_count"), 
                         col_types = "ccccc", show_col_types = FALSE)
  
  unique_names <- unique(refalt_char$name)
  cat("Unique names found:", length(unique_names), "\n")
  cat("First 10 names:", paste(head(unique_names, 10), collapse = ", "), "\n")
  
  # Check if we have the expected founders
  expected_founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
  found_founders <- intersect(expected_founders, unique_names)
  cat("Expected founders found:", paste(found_founders, collapse = ", "), "\n")
  cat("Missing founders:", paste(setdiff(expected_founders, unique_names), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("❌ Could not check sample names:", e$message, "\n")
})

cat("\n=== DIAGNOSIS COMPLETE ===\n")
