#!/usr/bin/env Rscript

# Centromere Wald Testing Analysis
# Command-line script to run complete centromere analysis

# Usage: Rscript run_centromere_analysis.R [results_file] [info_file] [output_file]
# Example: Rscript run_centromere_analysis.R combined_centromere_all_results.RDS info.ZINC2.txt Centromere_Wald_testing.RDS

library(tidyverse)

# Source the help functions
source("help_functions.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default values
results_file <- if (length(args) >= 1) args[1] else "combined_centromere_all_results.RDS"
info_file <- if (length(args) >= 2) args[2] else "info.ZINC2.txt"
output_file <- if (length(args) >= 3) args[3] else "Centromere_Wald_testing.RDS"

# Print usage if help requested
if (length(args) > 0 && (args[1] == "-h" || args[1] == "--help")) {
  cat("Usage: Rscript run_centromere_analysis.R [results_file] [info_file] [output_file]\n")
  cat("Default values:\n")
  cat("  results_file: combined_centromere_all_results.RDS\n")
  cat("  info_file: info.ZINC2.txt\n")
  cat("  output_file: Centromere_Wald_testing.RDS\n")
  quit(status = 0)
}

# Run the analysis
cat("Running centromere analysis with:\n")
cat("  Results file:", results_file, "\n")
cat("  Info file:", info_file, "\n")
cat("  Output file:", output_file, "\n\n")

tryCatch({
  results <- run_centromere_analysis(results_file, info_file, output_file)
  cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  quit(status = 1)
})
