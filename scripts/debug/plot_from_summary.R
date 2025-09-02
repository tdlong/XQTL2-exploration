#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript scripts/plot_from_summary.R <chr> <param_file> <output_dir> <sample_name> <center_pos_10kb> [width_10kb]")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
sample_name <- args[4]
center_pos_10kb <- as.numeric(args[5])
width_10kb <- if (length(args) >= 6) as.numeric(args[6]) else 20

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== PLOTTING FROM SUMMARY FILE ===\n")
cat("Chromosome:", chr, "\n")
cat("Sample:", sample_name, "\n")
cat("Center position (10kb):", center_pos_10kb, "\n")
cat("Width (10kb):", width_10kb, "\n\n")

# Load summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
if (!file.exists(summary_file)) {
  stop("Summary file not found. Run create_summary_file.R first: ", summary_file)
}

summary_data <- readRDS(summary_file)
cat("✓ Loaded summary file with", nrow(summary_data), "rows\n")

# Filter to sample and region
region_start_10kb <- center_pos_10kb - width_10kb/2
region_end_10kb <- center_pos_10kb + width_10kb/2

region_data <- summary_data %>%
  filter(sample == sample_name) %>%
  filter(pos_10kb >= region_start_10kb & pos_10kb <= region_end_10kb)

cat("✓ Filtered to", nrow(region_data), "rows for region", region_start_10kb, "-", region_end_10kb, "\n")

# Set color scheme: reds for fixed windows, greens for adaptive
method_colors <- c(
  "fixed_20kb" = "#8B0000",    # Dark red
  "fixed_50kb" = "#DC143C",    # Crimson
  "fixed_100kb" = "#FF6347",   # Tomato
  "fixed_200kb" = "#FF7F7F",   # Light coral
  "fixed_500kb" = "#FFB6C1",   # Light pink
  "adaptive_h4" = "#006400",   # Dark green
  "adaptive_h6" = "#228B22",   # Forest green
  "adaptive_h8" = "#32CD32",   # Lime green
  "adaptive_h10" = "#90EE90"   # Light green
)

# Create haplotype frequency plot (top panel)
p_haplo <- ggplot(region_data, aes(x = pos_10kb, y = B1, color = method)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = ","),
    breaks = seq(region_start_10kb, region_end_10kb, by = 2)  # Every 20kb
  ) +
  labs(
    title = paste("B1 Haplotype Frequencies (", width_10kb*10, "kb Zoom) -", sample_name),
    subtitle = paste("Chromosome:", chr, "Region:", format(region_start_10kb*10000, big.mark=","), "-", format(region_end_10kb*10000, big.mark=",")),
    x = NULL,  # No x-axis label for top panel
    y = "B1 Frequency",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none",  # Legend only on bottom panel
    axis.text.x = element_blank()  # No x-axis text for top panel
  )

# Create RMSE plot (bottom panel) - only for positions with SNPs
rmse_data <- region_data %>%
  filter(!is.na(rmse) & n_snps > 0)

if (nrow(rmse_data) > 0) {
  p_rmse <- ggplot(rmse_data, aes(x = pos_10kb, y = rmse, color = method)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(
      labels = function(x) format(x, scientific = FALSE, big.mark = ","),
      breaks = seq(region_start_10kb, region_end_10kb, by = 2)  # Every 20kb
    ) +
    labs(
      title = paste("SNP Imputation RMSE (", width_10kb*10, "kb Zoom) -", sample_name),
      subtitle = paste("Chromosome:", chr, "Region:", format(region_start_10kb*10000, big.mark=","), "-", format(region_end_10kb*10000, big.mark=","), "\nPositions with SNPs only"),
      x = "Position (10kb units)",
      y = "RMSE",
      color = "Method"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Create combined two-panel plot
  p_combined <- p_haplo / p_rmse + 
    plot_layout(heights = c(1, 1))
  
  # Save combined plot
  output_file <- file.path(results_dir, paste0("summary_plot_", chr, "_", sample_name, "_", center_pos_10kb, "_", width_10kb, ".png"))
  ggsave(output_file, p_combined, width = 14, height = 10, dpi = 300)
  cat("✓ Combined plot saved:", output_file, "\n")
  
} else {
  cat("⚠️  No positions with SNPs found in this region\n")
  # Save haplotype plot only
  output_file <- file.path(results_dir, paste0("summary_plot_", chr, "_", sample_name, "_", center_pos_10kb, "_", width_10kb, ".png"))
  ggsave(output_file, p_haplo, width = 12, height = 8, dpi = 300)
  cat("✓ Haplotype plot saved:", output_file, "\n")
}

# Print summary table for the region
cat("\n=== SUMMARY TABLE FOR REGION ===\n")
cat("Position (10kb) | Method        | B1     | Status | SNPs | RMSE\n")
cat("----------------|---------------|--------|--------|------|-----\n")

region_summary <- region_data %>%
  arrange(pos_10kb, method) %>%
  mutate(
    estimate_status = case_when(
      estimate_OK == 1 ~ "OK",
      estimate_OK == 0 ~ "FAIL",
      is.na(estimate_OK) ~ "NA"
    )
  )

for (i in 1:nrow(region_summary)) {
  row <- region_summary[i, ]
  pos_str <- sprintf("%-14.0f", row$pos_10kb)
  method_str <- sprintf("%-13s", row$method)
  b1_str <- sprintf("%-6.3f", row$B1)
  status_str <- sprintf("%-6s", row$estimate_status)
  snps_str <- sprintf("%-4d", row$n_snps)
  rmse_str <- ifelse(is.na(row$rmse), "NA", sprintf("%.3f", row$rmse))
  
  cat(pos_str, "|", method_str, "|", b1_str, "|", status_str, "|", snps_str, "|", rmse_str, "\n")
}

# Summary statistics
cat("\n=== REGION STATISTICS ===\n")
cat("Positions with SNPs:", sum(region_data$n_snps > 0), "/", nrow(region_data), "\n")
cat("Average SNPs per position (when > 0):", round(mean(region_data$n_snps[region_data$n_snps > 0], na.rm=TRUE), 1), "\n")

# Estimate_OK summary by method
ok_summary <- region_data %>%
  group_by(method) %>%
  summarise(
    total_positions = n(),
    ok_count = sum(estimate_OK == 1, na.rm = TRUE),
    fail_count = sum(estimate_OK == 0, na.rm = TRUE),
    na_count = sum(is.na(estimate_OK)),
    ok_rate = ok_count / total_positions,
    .groups = "drop"
  ) %>%
  arrange(desc(ok_rate))

cat("\nEstimate_OK summary by method:\n")
print(ok_summary)
