
#!/usr/bin/env Rscript

# =============================================================================
# Haplotype Method Evaluation Script - STREAMLINED VERSION
# =============================================================================
# 
# This script evaluates different haplotype estimation methods by comparing
# observed vs imputed SNP frequencies across different parameter combinations.
#
# METRICS CALCULATED:
# 1. Mean Squared Error (MSE) between observed and imputed frequencies
# 2. Coverage (proportion of positions with valid estimates)
# 3. Method comparison summary statistics
#
# USAGE:
# Rscript scripts/production/evaluate_imputation_methods.R <chr> <param_file> <output_dir>
#
# EXAMPLE:
# Rscript scripts/production/evaluate_imputation_methods.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript evaluate_imputation_methods.R <chr> <param_file> <output_dir>\n")
  cat("Example: Rscript evaluate_imputation_methods.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2\n")
  quit(status = 1)
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

cat("=== Haplotype Method Evaluation (Streamlined) ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameter file
source(param_file)

# Verify that names_in_bam is defined
if (!exists("names_in_bam")) {
  stop("Parameter file must define 'names_in_bam' variable")
}

cat("Samples to evaluate:", paste(names_in_bam, collapse = ", "), "\n")
cat("Total samples:", length(names_in_bam), "\n\n")

# Define euchromatin boundaries
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]

cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# =============================================================================
# Define parameter combinations to evaluate
# =============================================================================

fixed_params <- c(20, 50, 100, 200, 500)
adaptive_params <- c(4, 6, 8, 10)

estimators <- c(
  paste0("fixed_", fixed_params, "kb"),
  paste0("adaptive_h", adaptive_params)
)

cat("Evaluating", length(estimators), "estimators:\n")
cat("Fixed windows:", paste(fixed_params, "kb", collapse = ", "), "\n")
cat("Adaptive h_cutoffs:", paste(adaptive_params, collapse = ", "), "\n\n")

# =============================================================================
# Load observed SNP data from REFALT files
# =============================================================================

cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data not found: ", refalt_file)
}

# Load REFALT data
refalt_data <- read.table(refalt_file, header = TRUE)

# Transform to frequencies
cat("Converting counts to frequencies...\n")
refalt_processed <- refalt_data %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    total_count = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

# Filter for high-quality SNPs
cat("Filtering for high-quality SNPs...\n")
good_snps <- refalt_processed %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(total_count == 0),
    not_fixed = sum(total_count != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Get valid SNPs for evaluation (euchromatin only)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(observed_data %>% distinct(CHROM, POS)), "\n\n")

# Filter to euchromatin and only the samples we actually processed
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)  # Only process samples defined in parameter file

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# =============================================================================
# Function to load and process imputation results
# =============================================================================

load_imputation_results <- function(estimator, output_dir, chr) {
  results_dir <- file.path(output_dir, "haplotype_results")
  imputation_file <- file.path(results_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
  
  if (!file.exists(imputation_file)) {
    cat("  ⚠️  Imputation file not found:", imputation_file, "\n")
    return(NULL)
  }
  
  results <- read_rds(imputation_file)
  cat("  ✓ Loaded:", imputation_file, "(", nrow(results), "rows)\n")
  
  return(results)
}

# =============================================================================
# Load all imputation results
# =============================================================================

cat("Loading imputation results...\n")
imputation_results <- list()

for (estimator in estimators) {
  cat("Processing", estimator, "...\n")
  results <- load_imputation_results(estimator, output_dir, chr)
  if (!is.null(results)) {
    imputation_results[[estimator]] <- results
  }
}

cat("\n✓ Loaded", length(imputation_results), "imputation result files\n\n")

if (length(imputation_results) == 0) {
  stop("No imputation results found! Make sure to run SNP imputation first.")
}

# =============================================================================
# Calculate evaluation metrics
# =============================================================================

cat("Calculating evaluation metrics...\n")

evaluation_results <- list()

for (estimator in names(imputation_results)) {
  cat("Evaluating", estimator, "...\n")
  
  imputed_data <- imputation_results[[estimator]]
  
  # Merge observed and imputed data
  comparison_data <- observed_euchromatic %>%
    select(CHROM, POS, name, freq) %>%
    rename(observed_freq = freq) %>%
    left_join(
      imputed_data %>% 
        select(pos, sample, imputed_freq = imputed),
      by = c("POS" = "pos", "name" = "sample")
    )
  
  # Diagnostic information (simplified since we only process relevant samples)
  cat("  Sample coverage summary:\n")
  sample_coverage <- comparison_data %>%
    group_by(name) %>%
    summarize(
      n_total = n(),
      n_imputed = sum(!is.na(imputed_freq)),
      coverage_pct = round(n_imputed / n_total * 100, 1),
      .groups = "drop"
    )
  
  for (i in 1:nrow(sample_coverage)) {
    row <- sample_coverage[i, ]
    cat(sprintf("    %s: %d/%d SNPs (%.1f%%)\n", 
                row$name, row$n_imputed, row$n_total, row$coverage_pct))
  }
  
  # Calculate metrics (only for samples with imputed data)
  metrics <- comparison_data %>%
    group_by(name) %>%
    summarize(
      n_total = n(),
      n_imputed = sum(!is.na(imputed_freq)),
      coverage = n_imputed / n_total,
      .groups = "drop"
    ) %>%
    filter(n_imputed > 0) %>%  # Only include samples with some imputed data
    left_join(
      comparison_data %>% filter(!is.na(imputed_freq)),
      by = "name"
    ) %>%
    group_by(name, n_total, n_imputed, coverage) %>%
    reframe(
      mse = mean((observed_freq - imputed_freq)^2, na.rm = TRUE),
      rmse = sqrt(mse),
      mae = mean(abs(observed_freq - imputed_freq), na.rm = TRUE),
      correlation = ifelse(n_imputed > 1, 
                          tryCatch(cor(observed_freq, imputed_freq, use = "complete.obs"), 
                                  error = function(e) NA_real_), 
                          NA_real_)
    ) %>%
    mutate(estimator = estimator)
  
  evaluation_results[[estimator]] <- metrics
}

# Combine results
overall_metrics <- bind_rows(evaluation_results)

# =============================================================================
# Generate summary statistics
# =============================================================================

cat("\n=== OVERALL PERFORMANCE SUMMARY ===\n")

# Overall summary by estimator
overall_summary <- overall_metrics %>%
  group_by(estimator) %>%
  summarize(
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_mse = mean(mse, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_correlation = mean(correlation, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

# Create logical ordering: fixed windows (smallest to largest), then adaptive (smallest h to largest h)
estimator_order <- c(
  "fixed_20kb", "fixed_50kb", "fixed_100kb", "fixed_200kb", "fixed_500kb",
  "adaptive_h4", "adaptive_h6", "adaptive_h8", "adaptive_h10"
)

# Order by this logical sequence
overall_summary <- overall_summary %>%
  mutate(estimator = factor(estimator, levels = estimator_order)) %>%
  arrange(estimator)

# Format the table with better precision and column handling
cat("\nPerformance Summary (ordered logically: fixed windows 20kb→500kb, then adaptive h4→h10):\n")
cat("=", strrep("=", 100), "\n", sep = "")
cat(sprintf("%-20s %8s %10s %10s %10s %12s %8s\n", 
            "Estimator", "Coverage", "MSE", "RMSE", "MAE", "Correlation", "Samples"))
cat(strrep("-", 100), "\n")

for (i in 1:nrow(overall_summary)) {
  row <- overall_summary[i, ]
  cat(sprintf("%-20s %8.3f %10.6f %10.6f %10.6f %12.6f %8d\n",
              row$estimator,
              row$mean_coverage,
              row$mean_mse,
              row$mean_rmse,
              row$mean_mae,
              row$mean_correlation,
              row$n_samples))
}

cat(strrep("-", 100), "\n")
cat("\n")

# Add performance ranking information
cat("Performance Ranking (by MSE):\n")
performance_ranking <- overall_summary %>%
  arrange(mean_mse) %>%
  mutate(rank = row_number()) %>%
  select(rank, estimator, mean_mse, mean_correlation)

for (i in 1:nrow(performance_ranking)) {
  row <- performance_ranking[i, ]
  cat(sprintf("  %d. %s (MSE: %.6f, Correlation: %.3f)\n", 
              row$rank, row$estimator, row$mean_mse, row$mean_correlation))
}

cat("\n")

# =============================================================================
# Save results
# =============================================================================

cat("Saving results...\n")

# Create results subdirectory if it doesn't exist
results_dir <- file.path(output_dir, "haplotype_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Save summary statistics
summary_file <- file.path(results_dir, paste0("haplotype_evaluation_summary_", chr, ".RDS"))
saveRDS(overall_summary, summary_file)
cat("✓ Saved summary statistics:", summary_file, "\n")

# Save detailed metrics
detailed_file <- file.path(results_dir, paste0("haplotype_evaluation_detailed_", chr, ".RDS"))
saveRDS(overall_metrics, detailed_file)
cat("✓ Saved detailed metrics:", detailed_file, "\n")

# =============================================================================
# Print recommendations
# =============================================================================

cat("\n=== RECOMMENDATIONS ===\n")

# Best overall method (lowest MSE)
best_mse <- overall_summary %>% slice_min(mean_mse, n = 1)
cat("Best MSE performance:", best_mse$estimator, "(MSE =", round(best_mse$mean_mse, 6), ")\n")

# Best coverage
best_coverage <- overall_summary %>% slice_max(mean_coverage, n = 1)
cat("Best coverage:", best_coverage$estimator, "(coverage =", round(best_coverage$mean_coverage, 3), ")\n")

# Best correlation
best_corr <- overall_summary %>% slice_max(mean_correlation, n = 1)
cat("Best correlation:", best_corr$estimator, "(correlation =", round(best_corr$mean_correlation, 3), ")\n")

# Balanced performance (good coverage + low MSE)
balanced <- overall_summary %>%
  mutate(score = mean_correlation * mean_coverage / mean_mse) %>%
  slice_max(score, n = 1)
cat("Best balanced performance:", balanced$estimator, "(score =", round(balanced$score, 3), ")\n")

cat("\n=== OUTPUT FILES ===\n")
cat("Summary statistics:", summary_file, "\n")
cat("Detailed metrics:", detailed_file, "\n")

cat("\n=== EVALUATION COMPLETE ===\n")
cat("Note: This streamlined version focuses on core metrics only.\n")
cat("For advanced visualizations and sliding window analysis, use the next step in your pipeline.\n")
