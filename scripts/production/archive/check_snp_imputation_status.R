#!/usr/bin/env Rscript

# Check SNP Imputation Status
# Reads parameter table and checks which SNP imputation results exist

suppressPackageStartupMessages(library(tidyverse))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript check_snp_imputation_status.R <params_tsv> <output_dir>")
}

params_file <- args[1]
output_dir <- args[2]

cat("=== SNP IMPUTATION STATUS CHECK ===\n")
cat("Parameter file:", params_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Time:", format(Sys.time()), "\n\n")

# Read parameter combinations
if (!file.exists(params_file)) {
  stop("Parameter file not found: ", params_file)
}

params <- read.table(params_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(params) <- c("chromosome", "method", "parameter")

cat("=== PARAMETER COMBINATIONS ===\n")
cat("Total combinations:", nrow(params), "\n")
print(params)
cat("\n")

# Check which files exist
results_dir <- file.path(output_dir, "haplotype_results")
cat("=== CHECKING SNP IMPUTATION RESULTS ===\n")
cat("Results directory:", results_dir, "\n")
cat("Directory exists:", dir.exists(results_dir), "\n")
if (dir.exists(results_dir)) {
  cat("Files in directory:\n")
  files_in_dir <- list.files(results_dir, pattern = "snp_imputation.*\\.RDS$")
  if (length(files_in_dir) > 0) {
    for (f in files_in_dir) {
      cat("  ", f, "\n")
    }
  } else {
    cat("  No SNP imputation files found\n")
  }
}
cat("\n")

status_list <- list()

for (i in 1:nrow(params)) {
  chr <- params$chromosome[i]
  method <- params$method[i]
  param <- params$parameter[i]
  
  # Generate estimator name and file path
  if (method == "fixed") {
    estimator <- paste0("fixed_", param, "kb")
    expected_file <- file.path(results_dir, paste0("snp_imputation_fixed_", param, "kb_", chr, ".RDS"))
  } else if (method == "adaptive") {
    estimator <- paste0("adaptive_h", param)
    expected_file <- file.path(results_dir, paste0("snp_imputation_adaptive_h", param, "_", chr, ".RDS"))
  } else if (method == "smooth_h4") {
    estimator <- "smooth_h4"
    expected_file <- file.path(results_dir, paste0("snp_imputation_smooth_h4_", chr, ".RDS"))
  } else {
    stop("Unknown method: ", method)
  }
  
  # Check if file exists
  file_exists <- file.exists(expected_file)
  file_size <- if (file_exists) file.info(expected_file)$size else 0
  file_size_mb <- round(file_size / 1024^2, 2)
  
  # Store status
  status_list[[i]] <- list(
    combination = i,
    chromosome = chr,
    method = method,
    parameter = param,
    estimator = estimator,
    file_exists = file_exists,
    file_size_mb = file_size_mb,
    file_path = expected_file
  )
  
  # Print status
  status_symbol <- if (file_exists) "✓" else "❌"
  size_info <- if (file_exists) paste0(" (", file_size_mb, " MB)") else ""
  cat(sprintf("%s %s: %s%s\n", status_symbol, estimator, basename(expected_file), size_info))
}

cat("\n")

# Summary
completed <- sum(sapply(status_list, function(x) x$file_exists))
total <- length(status_list)
completion_rate <- round(completed / total * 100, 1)

cat("=== SUMMARY ===\n")
cat("Completed:", completed, "/", total, sprintf("(%.1f%%)\n", completion_rate))

if (completed < total) {
  missing <- total - completed
  cat("Missing:", missing, "files\n")
  
  cat("\nMissing combinations:\n")
  for (i in 1:length(status_list)) {
    if (!status_list[[i]]$file_exists) {
      cat("  ", status_list[[i]]$estimator, "\n")
    }
  }
  
  cat("\nTo run missing imputation jobs:\n")
  cat("sbatch scripts/production/snp_imputation_from_table.sh", params_file, "helpfiles/haplotype_parameters.R", output_dir, "\n")
} else {
  cat("🎉 All SNP imputation jobs complete!\n")
  
  total_size <- sum(sapply(status_list, function(x) x$file_size_mb))
  cat("Total data size:", round(total_size, 1), "MB\n")
  
  cat("\nNext steps:\n")
  cat("- Run evaluation: Rscript scripts/production/evaluate_imputation_methods.R\n")
  cat("- Check correlations and accuracy metrics\n")
}

cat("\n=== STATUS CHECK COMPLETE ===\n")
