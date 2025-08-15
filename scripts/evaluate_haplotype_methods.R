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
        select(pos, sample, imputed_freq = imputed),
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

# Format the table with better precision and column handling
cat("\nPerformance Summary (sorted by MSE):\n")
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
# Sliding Window Analysis
# =============================================================================

cat("\n=== SLIDING WINDOW ANALYSIS ===\n")

# Parameters for sliding window
window_size <- 250  # SNPs per window
step_size <- 50     # SNPs per step

# Get all unique SNP positions
all_positions <- observed_data %>%
  distinct(CHROM, POS) %>%
  arrange(POS) %>%
  pull(POS)

cat("Total SNP positions:", length(all_positions), "\n")
cat("Window size:", window_size, "SNPs\n")
cat("Step size:", step_size, "SNPs\n")

# Calculate window starts
window_starts <- seq(1, length(all_positions) - window_size + 1, by = step_size)
cat("Number of windows:", length(window_starts), "\n\n")

# Calculate sliding window performance for all estimators
cat("Calculating sliding windows for all estimators...\n")

# Create efficient sliding window analysis using vectorized operations
all_sliding_results <- list()

for (estimator in names(imputation_results)) {
  cat("Processing", estimator, "...\n")
  
  # Get imputed data for this estimator
  imputed_data <- imputation_results[[estimator]]
  
  # Create window assignments for all SNP positions
  window_assignments <- tibble(
    pos_idx = seq_along(all_positions),
    position = all_positions
  ) %>%
    mutate(
      window_start = ((pos_idx - 1) %/% step_size) * step_size + 1,
      window_id = (pos_idx - 1) %/% step_size + 1
    ) %>%
    filter(window_start <= length(all_positions) - window_size + 1) %>%
    mutate(
      window_end = window_start + window_size - 1,
      in_window = pos_idx >= window_start & pos_idx <= window_end
    ) %>%
    filter(in_window)
  
  # Get observed data with window assignments
  observed_with_windows <- observed_data %>%
    select(CHROM, POS, name, freq) %>%
    rename(observed_freq = freq) %>%
    left_join(
      window_assignments %>% select(position, window_id),
      by = c("POS" = "position")
    ) %>%
    filter(!is.na(window_id))
  
  # Get imputed data with window assignments
  imputed_with_windows <- imputed_data %>%
    select(pos, sample, imputed_freq = imputed) %>%
    left_join(
      window_assignments %>% select(position, window_id),
      by = c("pos" = "position")
    ) %>%
    filter(!is.na(window_id))
  
  # Calculate MSE for each window using vectorized operations
  window_performance <- observed_with_windows %>%
    left_join(
      imputed_with_windows,
      by = c("POS" = "pos", "name" = "sample", "window_id")
    ) %>%
    filter(!is.na(imputed_freq)) %>%
    group_by(window_id) %>%
    summarize(
      mse = mean((observed_freq - imputed_freq)^2, na.rm = TRUE),
      n_snps = n(),
      n_samples = n_distinct(name),
      .groups = "drop"
    ) %>%
    left_join(
      window_assignments %>%
        group_by(window_id) %>%
        summarize(
          mean_position = mean(position),
          .groups = "drop"
        ),
      by = "window_id"
    ) %>%
    mutate(estimator = estimator)
  
  all_sliding_results[[estimator]] <- window_performance
  
  cat("  ✓ Completed", estimator, "(", nrow(window_performance), "windows)\n")
}

# Combine all results
all_sliding_results <- bind_rows(all_sliding_results)

cat("\n✓ Sliding window analysis complete\n")
cat("Total windows analyzed:", nrow(all_sliding_results), "\n\n")

# Create sliding window plot
cat("Creating sliding window plot...\n")

sliding_plot <- all_sliding_results %>%
  ggplot(aes(x = mean_position / 1e6, y = sqrt(mse), color = estimator)) +
  geom_line(alpha = 0.8, linewidth = 0.8) +
  geom_point(alpha = 0.6, size = 0.5) +
  labs(
    title = paste("Sliding Window Performance:", chr),
    subtitle = paste("250 SNP windows, 50 SNP steps | Euchromatin:", euchromatin_start/1e6, "-", euchromatin_end/1e6, "Mb"),
    x = "Position (Mb)",
    y = "√MSE (Root Mean Squared Error)",
    color = "Estimator"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  ) +
  scale_color_viridis_d(option = "plasma")

# Create sliding window summary
cat("\n=== SLIDING WINDOW SUMMARY ===\n")

window_summary <- all_sliding_results %>%
  group_by(estimator) %>%
  summarize(
    mean_sqrt_mse = mean(sqrt(mse), na.rm = TRUE),
    median_sqrt_mse = median(sqrt(mse), na.rm = TRUE),
    sd_sqrt_mse = sd(sqrt(mse), na.rm = TRUE),
    min_sqrt_mse = min(sqrt(mse), na.rm = TRUE),
    max_sqrt_mse = max(sqrt(mse), na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_sqrt_mse)

# Print formatted sliding window summary
cat("\nSliding Window Summary (sorted by mean √MSE):\n")
cat("=", strrep("=", 120), "\n", sep = "")
cat(sprintf("%-20s %12s %12s %12s %12s %12s %8s\n", 
            "Estimator", "Mean √MSE", "Median √MSE", "SD √MSE", "Min √MSE", "Max √MSE", "Windows"))
cat(strrep("-", 120), "\n")

for (i in 1:nrow(window_summary)) {
  row <- window_summary[i, ]
  cat(sprintf("%-20s %12.6f %12.6f %12.6f %12.6f %12.6f %8d\n",
              row$estimator,
              row$mean_sqrt_mse,
              row$median_sqrt_mse,
              row$sd_sqrt_mse,
              row$min_sqrt_mse,
              row$max_sqrt_mse,
              row$n_windows))
}

cat(strrep("-", 120), "\n")
cat("\n")

# =============================================================================
# Save sliding window results (computationally expensive - save for future use)
# =============================================================================

cat("Saving sliding window results (computationally expensive analysis)...\n")

# Save sliding window results
sliding_file <- file.path(results_dir, paste0("sliding_window_results_", chr, ".RDS"))
saveRDS(all_sliding_results, sliding_file)
cat("✓ Saved sliding window data:", sliding_file, "\n")
cat("  (Use this file for future plot modifications without re-running analysis)\n")

# Save sliding window summary
sliding_summary_file <- file.path(results_dir, paste0("sliding_window_summary_", chr, ".RDS"))
saveRDS(window_summary, sliding_summary_file)
cat("✓ Saved sliding window summary:", sliding_summary_file, "\n")

# Save sliding window plot
sliding_plot_file <- file.path(results_dir, paste0("sliding_window_plot_", chr, ".png"))
ggsave(sliding_plot_file, sliding_plot, width = 14, height = 8, dpi = 300)
cat("✓ Saved sliding window plot:", sliding_plot_file, "\n")

cat("\n")

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
cat("Sliding window results:", sliding_file, "\n")
cat("Sliding window summary:", sliding_summary_file, "\n")
cat("Sliding window plot:", sliding_plot_file, "\n")

cat("\n=== EVALUATION COMPLETE ===\n")
