#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/production/plot_juice_vs_zinc2_comparison.R <chr> <juice_param_file> <zinc2_param_file> <method>")
}

chr <- args[1]
juice_param_file <- args[2]
zinc2_param_file <- args[3]
method_to_plot <- args[4]

cat("=== PLOTTING JUICE vs ZINC2 COMPARISON ===\n")
cat("Chromosome:", chr, "\n")
cat("JUICE parameter file:", juice_param_file, "\n")
cat("ZINC2 parameter file:", zinc2_param_file, "\n")
cat("Method:", method_to_plot, "\n\n")

# Load parameters for both datasets
source(juice_param_file)
juice_results_dir <- file.path("process/JUICE", "haplotype_results")

source(zinc2_param_file)
zinc2_results_dir <- file.path("process/ZINC2", "haplotype_results")

# Load summary files
juice_summary_file <- file.path(juice_results_dir, paste0("summary_", chr, ".RDS"))
zinc2_summary_file <- file.path(zinc2_results_dir, paste0("summary_", chr, ".RDS"))

if (!file.exists(juice_summary_file)) {
  stop("JUICE summary file not found: ", juice_summary_file)
}
if (!file.exists(zinc2_summary_file)) {
  stop("ZINC2 summary file not found: ", zinc2_summary_file)
}

cat("Loading summary files...\n")
juice_data <- readRDS(juice_summary_file)
zinc2_data <- readRDS(zinc2_summary_file)

cat("✓ JUICE data loaded:", nrow(juice_data), "rows\n")
cat("✓ ZINC2 data loaded:", nrow(zinc2_data), "rows\n")

# Get first sample from each dataset
juice_first_sample <- unique(juice_data$sample)[1]
zinc2_first_sample <- unique(zinc2_data$sample)[1]

cat("JUICE first sample:", juice_first_sample, "\n")
cat("ZINC2 first sample:", zinc2_first_sample, "\n\n")

# Filter data to specified method and first sample
juice_method_data <- juice_data %>%
  filter(sample == juice_first_sample) %>%
  filter(method == method_to_plot) %>%
  mutate(
    pos_10kb = pos / 10000,
    dataset = "JUICE"
  )

zinc2_method_data <- zinc2_data %>%
  filter(sample == zinc2_first_sample) %>%
  filter(method == method_to_plot) %>%
  mutate(
    pos_10kb = pos / 10000,
    dataset = "ZINC2"
  )

cat("JUICE", method_to_plot, "data points:", nrow(juice_method_data), "\n")
cat("ZINC2", method_to_plot, "data points:", nrow(zinc2_method_data), "\n\n")

if (nrow(juice_method_data) == 0) {
  stop("No JUICE data found for method: ", method_to_plot)
}
if (nrow(zinc2_method_data) == 0) {
  stop("No ZINC2 data found for method: ", method_to_plot)
}

# Combine datasets
combined_data <- bind_rows(juice_method_data, zinc2_method_data)

# Set color scheme
dataset_colors <- c(
  "JUICE" = "#4169E1",  # Royal blue
  "ZINC2" = "#DC143C"   # Crimson red
)

# Filter to only include reliable haplotype estimates
combined_reliable <- combined_data %>%
  filter(estimate_OK == TRUE)

cat("✓ Filtered to reliable estimates:\n")
cat("  JUICE:", sum(combined_reliable$dataset == "JUICE"), "positions\n")
cat("  ZINC2:", sum(combined_reliable$dataset == "ZINC2"), "positions\n\n")

# Create B1 frequency comparison plot (top panel)
# Set B1_freq to NA for unreliable estimates to create gaps in lines
combined_with_gaps <- combined_data %>%
  mutate(B1_freq = ifelse(estimate_OK == TRUE, B1_freq, NA))

p_b1 <- ggplot(combined_with_gaps, aes(x = pos_10kb, y = B1_freq, color = dataset)) +
  geom_point(size = 1, na.rm = TRUE, alpha = 0.7) +
  scale_color_manual(values = dataset_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(min(combined_data$pos_10kb)), round(max(combined_data$pos_10kb)), by = 100)  # Every 1Mb
  ) +
  labs(
    title = paste("B1 Haplotype Frequencies Comparison -", method_to_plot),
    subtitle = paste("Chromosome:", chr, "(positions in 10kb units)"),
    x = NULL,  # No x-axis label for top panel
    y = "B1 Frequency",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none",  # Legend only on bottom panel
    axis.text.x = element_blank()  # No x-axis text for top panel
  )

# Create MAE comparison plot (bottom panel) - smoothing curves only
p_mae <- ggplot(combined_data, aes(x = pos_10kb, y = MAE, color = dataset)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.05) +
  scale_color_manual(values = dataset_colors) +
  scale_y_continuous(
    breaks = seq(0, 0.30, by = 0.05),
    labels = sprintf("%.2f", seq(0, 0.30, by = 0.05)),
    limits = c(0, 0.30)
  ) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(min(combined_data$pos_10kb)), round(max(combined_data$pos_10kb)), by = 100)  # Every 1Mb
  ) +
  labs(
    title = paste("SNP Imputation Mean Absolute Error Comparison -", method_to_plot),
    x = "Position (10kb units)",
    y = "Mean Absolute Error",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Create combined two-panel plot
p_combined <- p_b1 / p_mae + 
  plot_layout(heights = c(1, 1.2))  # Slightly taller bottom panel for legend

# Save combined plot
plot_file <- file.path(juice_results_dir, paste0("juice_vs_zinc2_comparison_", chr, "_", method_to_plot, ".png"))
ggsave(plot_file, p_combined, width = 16, height = 10, dpi = 300)

cat("✓ Combined comparison plot saved to:", plot_file, "\n")

# Print summary statistics
cat("\n=== COMPARISON SUMMARY ===\n")

# B1 frequency comparison
b1_summary <- combined_reliable %>%
  group_by(dataset) %>%
  summarize(
    mean_b1 = mean(B1_freq, na.rm = TRUE),
    median_b1 = median(B1_freq, na.rm = TRUE),
    min_b1 = min(B1_freq, na.rm = TRUE),
    max_b1 = max(B1_freq, na.rm = TRUE),
    .groups = "drop"
  )

cat("B1 Frequency Statistics:\n")
for (i in 1:nrow(b1_summary)) {
  row <- b1_summary[i, ]
  cat(sprintf("  %s: mean=%.3f, median=%.3f, range=[%.3f-%.3f]\n",
              row$dataset, row$mean_b1, row$median_b1, row$min_b1, row$max_b1))
}

# MAE comparison
mae_summary <- combined_data %>%
  group_by(dataset) %>%
  summarize(
    mean_mae = mean(MAE, na.rm = TRUE),
    median_mae = median(MAE, na.rm = TRUE),
    min_mae = min(MAE, na.rm = TRUE),
    max_mae = max(MAE, na.rm = TRUE),
    positions_with_mae = sum(!is.na(MAE)),
    total_positions = n(),
    mae_coverage_pct = round(positions_with_mae / total_positions * 100, 1),
    .groups = "drop"
  )

cat("\nMAE Statistics:\n")
for (i in 1:nrow(mae_summary)) {
  row <- mae_summary[i, ]
  cat(sprintf("  %s: mean=%.3f, median=%.3f, range=[%.3f-%.3f], coverage=%d/%d (%.1f%%)\n",
              row$dataset, row$mean_mae, row$median_mae, row$min_mae, row$max_mae,
              row$positions_with_mae, row$total_positions, row$mae_coverage_pct))
}

cat("\n✓ Comparison plotting complete!\n")
