#!/usr/bin/env Rscript

# =============================================================================
# Haplotype Method Evaluation Script
# =============================================================================
# 
# This script evaluates different haplotype estimation methods by comparing
# observed vs imputed SNP frequencies across different parameter combinations.
#
# METRICS CALCULATED:
# 1. Mean Squared Error (MSE) between observed and imputed frequencies
# 2. Coverage (proportion of positions with valid estimates)
# 3. Regional performance (MSE by chromosomal windows)
# 4. Method comparison summary statistics
#
# USAGE:
# Rscript scripts/evaluate_haplotype_methods.R <chr> <param_file> <output_dir> [window_size_kb]
#
# EXAMPLE:
# Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/haplotype_parameters.R process/JUICE 1000
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(patchwork)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript evaluate_haplotype_methods.R <chr> <param_file> <output_dir> [window_size_kb]\n")
  cat("Example: Rscript evaluate_haplotype_methods.R chr2R helpfiles/haplotype_parameters.R process/JUICE 1000\n")
  quit(status = 1)
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
window_size_kb <- ifelse(length(args) >= 4, as.numeric(args[4]), 1000)

cat("=== Haplotype Method Evaluation ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Regional window size:", window_size_kb, "kb\n\n")

# Load parameter file
source(param_file)

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

# Load REFALT data - use same format as other scripts
refalt_data <- read.table(refalt_file, header = TRUE)

# Transform to frequencies (same as other scripts)
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

# Filter for high-quality SNPs (same as other scripts)
cat("Filtering for high-quality SNPs...\n")
good_snps <- refalt_processed %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(total_count == 0),
    not_fixed = sum(total_count != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Get valid SNPs for evaluation (euchromatin only)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(observed_data %>% distinct(CHROM, POS)), "\n\n")

# Filter to euchromatin and non-founder samples
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(!name %in% founders)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# =============================================================================
# Function to load and process imputation results
# =============================================================================

load_imputation_results <- function(estimator, output_dir, chr) {
  # Look in the haplotype_results subdirectory
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
        select(pos, sample, imputed_freq = freq),
      by = c("POS" = "pos", "name" = "sample")
    )
  
  # Calculate metrics
  metrics <- comparison_data %>%
    group_by(name) %>%
    summarize(
      n_total = n(),
      n_imputed = sum(!is.na(imputed_freq)),
      coverage = n_imputed / n_total,
      mse = mean((observed_freq - imputed_freq)^2, na.rm = TRUE),
      rmse = sqrt(mse),
      mae = mean(abs(observed_freq - imputed_freq), na.rm = TRUE),
      correlation = cor(observed_freq, imputed_freq, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(estimator = estimator)
  
  evaluation_results[[estimator]] <- metrics
  
  # Add regional analysis
  regional_data <- comparison_data %>%
    filter(!is.na(imputed_freq)) %>%
    mutate(
      window_start = floor(POS / (window_size_kb * 1000)) * (window_size_kb * 1000),
      window_end = window_start + (window_size_kb * 1000),
      window_mid = (window_start + window_end) / 2
    ) %>%
    group_by(window_mid) %>%
    summarize(
      n_snps = n(),
      regional_mse = mean((observed_freq - imputed_freq)^2, na.rm = TRUE),
      regional_rmse = sqrt(regional_mse),
      regional_mae = mean(abs(observed_freq - imputed_freq), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(estimator = estimator)
  
  evaluation_results[[paste0(estimator, "_regional")]] <- regional_data
}

# Combine results
overall_metrics <- bind_rows(evaluation_results[!grepl("_regional$", names(evaluation_results))])
regional_metrics <- bind_rows(evaluation_results[grepl("_regional$", names(evaluation_results))])

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
  ) %>%
  arrange(mean_mse)

print(overall_summary)

# =============================================================================
# Create visualizations
# =============================================================================

cat("\nGenerating visualizations...\n")

# 1. MSE comparison across methods
p1 <- overall_metrics %>%
  ggplot(aes(x = reorder(estimator, mse), y = mse, fill = estimator)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Mean Squared Error by Haplotype Estimation Method",
    x = "Estimator",
    y = "MSE (Observed vs Imputed Frequencies)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Coverage comparison
p2 <- overall_metrics %>%
  ggplot(aes(x = reorder(estimator, coverage), y = coverage, fill = estimator)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Coverage by Haplotype Estimation Method",
    x = "Estimator", 
    y = "Proportion of Positions with Estimates"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Regional performance (MSE across chromosome)
p3 <- regional_metrics %>%
  ggplot(aes(x = window_mid / 1e6, y = regional_rmse, color = estimator)) +
  geom_line(alpha = 0.7) +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(
    title = "Regional Performance: RMSE Across Chromosome",
    x = "Position (Mb)",
    y = "RMSE"
  ) +
  theme_minimal() +
  facet_wrap(~estimator, scales = "free_y", ncol = 3)

# 4. Correlation vs Coverage trade-off
p4 <- overall_metrics %>%
  ggplot(aes(x = coverage, y = correlation, color = estimator)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.5, size = 2.5) +
  labs(
    title = "Coverage vs Correlation Trade-off",
    x = "Coverage",
    y = "Correlation with Observed Frequencies"
  ) +
  theme_minimal()

# Combine plots
combined_plot <- (p1 + p2) / (p3) / (p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste("Haplotype Method Evaluation:", chr),
    subtitle = paste("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp")
  )

# =============================================================================
# Save results
# =============================================================================

# Create results subdirectory if it doesn't exist
results_dir <- file.path(output_dir, "haplotype_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Save summary statistics
summary_file <- file.path(results_dir, paste0("haplotype_evaluation_summary_", chr, ".RDS"))
saveRDS(overall_summary, summary_file)

# Save detailed metrics
detailed_file <- file.path(results_dir, paste0("haplotype_evaluation_detailed_", chr, ".RDS"))
saveRDS(overall_metrics, detailed_file)

# Save regional analysis
regional_file <- file.path(results_dir, paste0("haplotype_evaluation_regional_", chr, ".RDS"))
saveRDS(regional_metrics, regional_file)

# Save plot
plot_file <- file.path(results_dir, paste0("haplotype_evaluation_plots_", chr, ".png"))
ggsave(plot_file, combined_plot, width = 16, height = 12, dpi = 300)

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
cat("Regional analysis:", regional_file, "\n")
cat("Evaluation plots:", plot_file, "\n")

cat("\n=== EVALUATION COMPLETE ===\n")
