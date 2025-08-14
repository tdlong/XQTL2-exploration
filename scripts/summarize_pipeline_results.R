#!/usr/bin/env Rscript

# =============================================================================
# Pipeline Results Summary
# =============================================================================
# 
# This script summarizes the results of the XQTL2 pipeline by:
# 1. Reading all haplotype estimation files and calculating success rates
# 2. Reading all SNP imputation files and calculating success rates
# 3. Reporting results by parameter combination
#
# USAGE:
# Rscript scripts/summarize_pipeline_results.R <output_dir> <chromosome>
#
# EXAMPLE:
# Rscript scripts/summarize_pipeline_results.R process/JUICE chr2R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(glue)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript scripts/summarize_pipeline_results.R <output_dir> <chromosome>\n")
  cat("Example: Rscript scripts/summarize_pipeline_results.R process/JUICE chr2R\n")
  quit(status = 1)
}

output_dir <- args[1]
chr <- args[2]

cat("=== PIPELINE RESULTS SUMMARY ===\n")
cat("Output directory:", output_dir, "\n")
cat("Chromosome:", chr, "\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# Load parameter combinations
# =============================================================================

cat("=== LOADING PARAMETER COMBINATIONS ===\n")

# Read the parameter file to know what combinations to expect
# Extract chromosome number (e.g., chr2R -> 2R)
chr_num <- gsub("chr", "", chr)
param_file <- file.path("helpfiles", glue("haplotype_params.{chr_num}.tsv"))
if (!file.exists(param_file)) {
  cat("❌ Parameter file not found:", param_file, "\n")
  quit(status = 1)
}

params <- read_tsv(param_file, col_names = c("chromosome", "method", "parameter"), col_types = "ccc")
cat("✓ Parameter file loaded:", nrow(params), "combinations\n")

# Create expected file patterns
expected_files <- params %>%
  mutate(
    file_pattern = case_when(
      method == "fixed" ~ glue("fixed_window_{parameter}kb_results_{chr}.RDS"),
      method == "adaptive" ~ glue("adaptive_window_h{parameter}_results_{chr}.RDS")
    ),
    estimator = case_when(
      method == "fixed" ~ glue("fixed_{parameter}kb"),
      method == "adaptive" ~ glue("adaptive_h{parameter}")
    )
  )

cat("Expected haplotype files:\n")
for (i in 1:nrow(expected_files)) {
  cat("  ", expected_files$file_pattern[i], "\n")
}

# =============================================================================
# Analyze haplotype estimation results
# =============================================================================

cat("\n=== HAPLOTYPE ESTIMATION RESULTS ===\n")

haplotype_dir <- file.path(output_dir, "haplotype_results")
if (!dir.exists(haplotype_dir)) {
  cat("❌ Haplotype results directory not found:", haplotype_dir, "\n")
  quit(status = 1)
}

haplotype_summary <- list()

for (i in 1:nrow(expected_files)) {
  file_path <- file.path(haplotype_dir, expected_files$file_pattern[i])
  method <- expected_files$method[i]
  parameter <- expected_files$parameter[i]
  estimator <- expected_files$estimator[i]
  
  if (file.exists(file_path)) {
    cat("Analyzing:", expected_files$file_pattern[i], "\n")
    
    # Load the file
    results <- readRDS(file_path)
    
    # Identify founder columns (exclude chr, pos, sample)
    non_founder_cols <- c("chr", "pos", "sample")
    founder_cols <- names(results)[!names(results) %in% non_founder_cols]
    
    cat("  Founder columns found:", paste(founder_cols, collapse = ", "), "\n")
    
    # Calculate success rate (rows where at least one founder has non-NA value)
    total_rows <- nrow(results)
    successful_rows <- sum(rowSums(!is.na(results[founder_cols])) > 0)
    success_rate <- successful_rows / total_rows * 100
    
    # Calculate sample-level statistics
    sample_stats <- results %>%
      group_by(sample) %>%
      summarize(
        total_positions = n(),
        successful_positions = sum(rowSums(!is.na(select(., all_of(founder_cols)))) > 0),
        success_rate = successful_positions / total_positions * 100
      )
    
    haplotype_summary[[estimator]] <- list(
      method = method,
      parameter = parameter,
      file_path = file_path,
      total_rows = total_rows,
      successful_rows = successful_rows,
      success_rate = success_rate,
      founder_cols = founder_cols,
      sample_stats = sample_stats,
      file_size_mb = round(file_info(file_path)$size / 1024^2, 1)
    )
    
    cat("  Success rate:", round(success_rate, 1), "% (", successful_rows, "/", total_rows, ")\n")
    cat("  File size:", haplotype_summary[[estimator]]$file_size_mb, "MB\n")
    
  } else {
    cat("❌ File not found:", expected_files$file_pattern[i], "\n")
  }
}

# =============================================================================
# Analyze SNP imputation results
# =============================================================================

cat("\n=== SNP IMPUTATION RESULTS ===\n")

snp_imputation_summary <- list()

for (i in 1:nrow(expected_files)) {
  estimator <- expected_files$estimator[i]
  file_pattern <- glue("snp_imputation_{estimator}_{chr}.RDS")
  file_path <- file.path(haplotype_dir, file_pattern)
  
  if (file.exists(file_path)) {
    cat("Analyzing:", file_pattern, "\n")
    
    # Load the file
    results <- readRDS(file_path)
    
    # Calculate success rate (rows where imputed is not NA)
    total_rows <- nrow(results)
    successful_rows <- sum(!is.na(results$imputed))
    success_rate <- successful_rows / total_rows * 100
    
    # Calculate sample-level statistics
    sample_stats <- results %>%
      group_by(sample) %>%
      summarize(
        total_positions = n(),
        successful_positions = sum(!is.na(imputed)),
        success_rate = successful_positions / total_positions * 100,
        mean_observed = mean(observed, na.rm = TRUE),
        mean_imputed = mean(imputed, na.rm = TRUE),
        correlation = cor(observed, imputed, use = "complete.obs")
      )
    
    snp_imputation_summary[[estimator]] <- list(
      file_path = file_path,
      total_rows = total_rows,
      successful_rows = successful_rows,
      success_rate = success_rate,
      sample_stats = sample_stats,
      file_size_mb = round(file_info(file_path)$size / 1024^2, 1)
    )
    
    cat("  Success rate:", round(success_rate, 1), "% (", successful_rows, "/", total_rows, ")\n")
    cat("  File size:", snp_imputation_summary[[estimator]]$file_size_mb, "MB\n")
    
  } else {
    cat("❌ File not found:", file_pattern, "\n")
  }
}

# =============================================================================
# Generate summary report
# =============================================================================

cat("\n=== SUMMARY REPORT ===\n")

# Haplotype estimation summary
if (length(haplotype_summary) > 0) {
  cat("\nHaplotype Estimation Results:\n")
  cat("=", strrep("=", 50), "\n")
  
  # Create summary table
  haplotype_table <- data.frame(
    Method = sapply(haplotype_summary, function(x) x$method),
    Parameter = sapply(haplotype_summary, function(x) x$parameter),
    Success_Rate = sapply(haplotype_summary, function(x) round(x$success_rate, 1)),
    Total_Rows = sapply(haplotype_summary, function(x) x$total_rows),
    Successful_Rows = sapply(haplotype_summary, function(x) x$successful_rows),
    File_Size_MB = sapply(haplotype_summary, function(x) x$file_size_mb)
  )
  
  print(haplotype_table)
  
  # Overall statistics
  total_haplotype_files <- length(haplotype_summary)
  avg_success_rate <- mean(sapply(haplotype_summary, function(x) x$success_rate))
  cat("\nOverall haplotype estimation:\n")
  cat("  Files completed:", total_haplotype_files, "/", nrow(expected_files), "\n")
  cat("  Average success rate:", round(avg_success_rate, 1), "%\n")
  
  # Method comparison
  if (length(unique(haplotype_table$Method)) > 1) {
    cat("\nMethod comparison:\n")
    method_comparison <- haplotype_table %>%
      group_by(Method) %>%
      summarize(
        avg_success_rate = mean(Success_Rate),
        n_files = n()
      )
    print(method_comparison)
  }
}

# SNP imputation summary
if (length(snp_imputation_summary) > 0) {
  cat("\nSNP Imputation Results:\n")
  cat("=", strrep("=", 50), "\n")
  
  # Create summary table
  snp_table <- data.frame(
    Estimator = names(snp_imputation_summary),
    Success_Rate = sapply(snp_imputation_summary, function(x) round(x$success_rate, 1)),
    Total_Rows = sapply(snp_imputation_summary, function(x) x$total_rows),
    Successful_Rows = sapply(snp_imputation_summary, function(x) x$successful_rows),
    File_Size_MB = sapply(snp_imputation_summary, function(x) x$file_size_mb)
  )
  
  print(snp_table)
  
  # Overall statistics
  total_snp_files <- length(snp_imputation_summary)
  avg_snp_success_rate <- mean(sapply(snp_imputation_summary, function(x) x$success_rate))
  cat("\nOverall SNP imputation:\n")
  cat("  Files completed:", total_snp_files, "/", nrow(expected_files), "\n")
  cat("  Average success rate:", round(avg_snp_success_rate, 1), "%\n")
}

# =============================================================================
# Recommendations
# =============================================================================

cat("\n=== RECOMMENDATIONS ===\n")

if (length(haplotype_summary) >= 2) {
  cat("✅ Ready for haplotype method comparison:\n")
  param_file <- "helpfiles/JUICE/JUICE_haplotype_parameters.R"
  cat("  Rscript scripts/evaluate_haplotype_methods.R", chr, param_file, output_dir, "\n")
}

if (length(haplotype_summary) > 0) {
  cat("\n✅ Individual file validation:\n")
  for (estimator in names(haplotype_summary)) {
    cat("  Rscript scripts/peek_haplotype_results.R", haplotype_summary[[estimator]]$file_path, "\n")
  }
}

cat("\n=== DONE ===\n")
