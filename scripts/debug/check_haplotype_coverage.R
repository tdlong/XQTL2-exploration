#!/usr/bin/env Rscript

# Check Haplotype Coverage
# Summarizes which haplotype files exist and their estimate_OK success rates

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("=== HAPLOTYPE COVERAGE SUMMARY ===\n\n")

# Parameters
chr <- "chr2R"
results_dir <- "process/JUICE/haplotype_results"

# Expected haplotype file patterns
expected_files <- data.frame(
  method = c(rep("fixed", 5), rep("adaptive", 4)),
  parameter = c("20kb", "50kb", "100kb", "200kb", "500kb", "h4", "h6", "h8", "h10"),
  file_pattern = c(
    paste0("fixed_window_20kb_results_", chr, ".RDS"),
    paste0("fixed_window_50kb_results_", chr, ".RDS"),
    paste0("fixed_window_100kb_results_", chr, ".RDS"),
    paste0("fixed_window_200kb_results_", chr, ".RDS"),
    paste0("fixed_window_500kb_results_", chr, ".RDS"),
    paste0("adaptive_window_h4_results_", chr, ".RDS"),
    paste0("adaptive_window_h6_results_", chr, ".RDS"),
    paste0("adaptive_window_h8_results_", chr, ".RDS"),
    paste0("adaptive_window_h10_results_", chr, ".RDS")
  ),
  stringsAsFactors = FALSE
)

cat("Checking for haplotype result files in:", results_dir, "\n\n")

# Check which files exist
file_summary <- expected_files %>%
  mutate(
    file_path = file.path(results_dir, file_pattern),
    file_exists = file.exists(file_path),
    estimator = paste0(method, "_", parameter)
  )

# Show file existence summary
cat("=== FILE EXISTENCE ===\n")
for (i in 1:nrow(file_summary)) {
  status <- ifelse(file_summary$file_exists[i], "✓", "✗")
  cat(sprintf("%-15s %s %s\n", file_summary$estimator[i], status, file_summary$file_pattern[i]))
}

existing_files <- file_summary %>% filter(file_exists)
cat("\nFound", nrow(existing_files), "of", nrow(expected_files), "expected files\n\n")

if (nrow(existing_files) == 0) {
  cat("❌ No haplotype result files found!\n")
  quit(status = 1)
}

# Analyze each existing file
cat("=== HAPLOTYPE COVERAGE ANALYSIS ===\n")

coverage_summary <- data.frame()

for (i in 1:nrow(existing_files)) {
  estimator <- existing_files$estimator[i]
  file_path <- existing_files$file_path[i]
  
  cat("\nAnalyzing:", estimator, "\n")
  
  # Load the results
  results <- readRDS(file_path)
  
  # Show basic structure
  cat("  Data structure:", nrow(results), "rows,", ncol(results), "columns\n")
  cat("  Columns:", paste(names(results), collapse = ", "), "\n")
  
  # Check for estimate_OK column
  if ("estimate_OK" %in% names(results)) {
    
    # Overall summary
    total_positions <- nrow(results)
    valid_estimates <- sum(results$estimate_OK == 1, na.rm = TRUE)
    failed_estimates <- sum(results$estimate_OK == 0, na.rm = TRUE)
    na_estimates <- sum(is.na(results$estimate_OK))
    
    cat("  Total positions:", total_positions, "\n")
    cat("  Valid estimates (estimate_OK = 1):", valid_estimates, sprintf("(%.1f%%)", valid_estimates/total_positions*100), "\n")
    cat("  Failed estimates (estimate_OK = 0):", failed_estimates, sprintf("(%.1f%%)", failed_estimates/total_positions*100), "\n")
    cat("  NA estimates:", na_estimates, sprintf("(%.1f%%)", na_estimates/total_positions*100), "\n")
    
    # Per-sample summary if sample column exists
    if ("sample" %in% names(results)) {
      sample_summary <- results %>%
        group_by(sample) %>%
        summarise(
          total_pos = n(),
          valid_pos = sum(estimate_OK == 1, na.rm = TRUE),
          failed_pos = sum(estimate_OK == 0, na.rm = TRUE),
          na_pos = sum(is.na(estimate_OK)),
          success_rate = round(valid_pos / total_pos * 100, 1),
          .groups = "drop"
        )
      
      cat("  Per-sample coverage:\n")
      for (j in 1:nrow(sample_summary)) {
        cat(sprintf("    %-8s: %4d/%4d valid (%.1f%%)\n", 
                   sample_summary$sample[j], 
                   sample_summary$valid_pos[j], 
                   sample_summary$total_pos[j],
                   sample_summary$success_rate[j]))
      }
    }
    
    # Store summary for comparison
    if ("sample" %in% names(results)) {
      # Multi-sample file
      sample_data <- results %>%
        group_by(sample) %>%
        summarise(
          total_positions = n(),
          valid_estimates = sum(estimate_OK == 1, na.rm = TRUE),
          success_rate = round(valid_estimates / total_positions * 100, 1),
          .groups = "drop"
        ) %>%
        mutate(estimator = estimator, method = existing_files$method[i], parameter = existing_files$parameter[i])
      
      coverage_summary <- bind_rows(coverage_summary, sample_data)
    } else {
      # Single summary
      single_data <- data.frame(
        estimator = estimator,
        method = existing_files$method[i],
        parameter = existing_files$parameter[i],
        sample = "combined",
        total_positions = total_positions,
        valid_estimates = valid_estimates,
        success_rate = round(valid_estimates / total_positions * 100, 1)
      )
      coverage_summary <- bind_rows(coverage_summary, single_data)
    }
    
  } else {
    cat("  ⚠️  No 'estimate_OK' column found - this may be old format\n")
  }
}

# Final comparison table
if (nrow(coverage_summary) > 0) {
  cat("\n=== FINAL SUMMARY TABLE ===\n")
  
  # Pivot wider for easier comparison
  summary_table <- coverage_summary %>%
    select(estimator, sample, success_rate) %>%
    pivot_wider(names_from = sample, values_from = success_rate, values_fill = 0)
  
  print(summary_table)
  
  cat("\n=== METHOD COMPARISON ===\n")
  method_comparison <- coverage_summary %>%
    group_by(estimator, method, parameter) %>%
    summarise(
      avg_success_rate = round(mean(success_rate), 1),
      total_samples = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_success_rate))
  
  print(method_comparison)
}

cat("\n=== HAPLOTYPE COVERAGE CHECK COMPLETE ===\n")
