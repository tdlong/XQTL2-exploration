#!/usr/bin/env Rscript

# Plot summary for entire chromosome - adaptive_h4 only
# Usage: Rscript scripts/production/plot_summary_chromosome.R <chr> <params_file> <results_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/production/plot_summary_chromosome.R <chr> <params_file> <results_dir>")
}

chr <- args[1]
params_file <- args[2]
results_dir <- args[3]

cat("=== PLOTTING SUMMARY CHROMOSOME ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", results_dir, "\n\n")

# Load haplotype parameters
if (!file.exists(params_file)) {
  stop("Parameters file not found:", params_file)
}

source(params_file)

# Check if results directory exists
if (!dir.exists(results_dir)) {
  stop("Results directory not found:", results_dir)
}

# Load summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
if (!file.exists(summary_file)) {
  stop("Summary file not found:", summary_file)
}

cat("Loading summary file...\n")
summary_data <- readRDS(summary_file)
cat("✓ Summary data loaded:", nrow(summary_data), "rows\n")

# Get the first sample
first_sample <- unique(summary_data$sample)[1]
cat("Using first sample:", first_sample, "\n\n")

# Filter data to only adaptive_h4 and first sample (no region filtering - entire chromosome)
region_data <- summary_data %>%
  filter(sample == first_sample) %>%
  filter(method == "adaptive_h4") %>%
  mutate(pos_10kb = pos / 10000)  # Convert to 10kb units for plotting

cat("Data points for adaptive_h4:", nrow(region_data), "\n")

if (nrow(region_data) == 0) {
  stop("No adaptive_h4 data found")
}

# Show position range in data
pos_range <- range(region_data$pos)
cat("Position range in data:", format(pos_range[1], big.mark=","), "to", format(pos_range[2], big.mark=","), "base pairs\n")
cat("Position range in 10kb units:", round(pos_range[1]/10000), "to", round(pos_range[2]/10000), "\n\n")

# Set color scheme for adaptive_h4 only
method_colors <- c("adaptive_h4" = "#006400")  # Dark green

# Filter data to only include reliable haplotype estimates
region_data_reliable <- region_data %>%
  filter(estimate_OK == TRUE)

cat("✓ Filtered to reliable estimates:", nrow(region_data_reliable), "rows\n")
cat("  Original data:", nrow(region_data), "rows\n")
cat("  Reliable estimates:", nrow(region_data_reliable), "rows\n")
cat("  Unreliable estimates:", nrow(region_data) - nrow(region_data_reliable), "rows\n")

# Create haplotype frequency plot (top panel)
# Set B1_freq to NA for unreliable estimates to create gaps in lines
region_data_with_gaps <- region_data %>%
  mutate(B1_freq = ifelse(estimate_OK == TRUE, B1_freq, NA))

p_haplo <- ggplot(region_data_with_gaps, aes(x = pos_10kb, y = B1_freq, color = method)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(pos_range[1]/10000), round(pos_range[2]/10000), by = 100)  # Every 1Mb (100 * 10kb)
  ) +
  labs(
    title = paste("B1 Haplotype Frequencies -", first_sample),
    subtitle = paste("Chromosome:", chr, "(positions in 10kb units)"),
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

# Create MAE plot (middle panel) - show all positions, NA means no imputation possible
p_mae <- ggplot(region_data, aes(x = pos_10kb, y = MAE, color = method)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = method_colors) +
  scale_y_continuous(
    breaks = seq(0.07, 0.15, by = 0.01),
    labels = sprintf("%.2f", seq(0.07, 0.15, by = 0.01)),
    limits = c(0.07, 0.15)
  ) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(pos_range[1]/10000), round(pos_range[2]/10000), by = 100)  # Every 1Mb (100 * 10kb)
  ) +
  labs(
    title = paste("SNP Imputation Mean Absolute Error -", first_sample),
    x = NULL,  # No x-axis label for middle panel
    y = "Mean Absolute Error",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none",  # Legend only on bottom panel
    axis.text.x = element_blank()  # No x-axis text for middle panel
  )

# Create SNP count plot (bottom panel) - identical structure to panel 1, just different y-axis
p_snps <- ggplot(region_data, aes(x = pos_10kb, y = NSNPs, color = method)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(pos_range[1]/10000), round(pos_range[2]/10000), by = 100)  # Every 1Mb (100 * 10kb)
  ) +
  # Use log scale for y-axis with reasonable number of labels
  scale_y_log10(
    breaks = c(100, 200, 500, 1000, 2000, 5000),
    labels = c("100", "200", "500", "1000", "2000", "5000")
  ) +
  labs(
    title = paste("Number of SNPs per Window -", first_sample),
    subtitle = paste("Chromosome:", chr, "(SNPs actually used in haplotype estimation: from haplotype files)"),
    x = "Position (10kb units)",
    y = "Number of SNPs (log scale)",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Create combined three-panel plot
p_combined <- p_haplo / p_mae / p_snps + 
  plot_layout(heights = c(1, 1, 1.2))  # Slightly taller bottom panel for legend

# Save combined plot
plot_file <- file.path(results_dir, "haplotype_results", paste0("chromosome_summary_", chr, "_", first_sample, "_adaptive_h4.png"))
ggsave(plot_file, p_combined, width = 16, height = 12, dpi = 300)

cat("✓ Combined chromosome plot saved to:", plot_file, "\n")

# Print summary table for the chromosome
cat("\n=== SUMMARY TABLE FOR CHROMOSOME ===\n")
cat("Pos (10kb) | Pos (bp)      | B1     | Status | SNPs\n")
cat("------------|---------------|--------|--------|------\n")

# Show first 10 and last 10 positions for overview
first_10 <- head(region_data, 10)
last_10 <- tail(region_data, 10)

for (i in 1:nrow(first_10)) {
  row <- first_10[i, ]
  pos_10kb_str <- sprintf("%-10.0f", row$pos_10kb)
  pos_bp_str <- sprintf("%-13.0f", row$pos)
  b1_str <- sprintf("%-6.3f", row$B1_freq)
  status_str <- sprintf("%-6s", ifelse(row$estimate_OK == 1, "OK", "FAIL"))
  snps_str <- sprintf("%-4d", row$NSNPs)
  
  cat(pos_10kb_str, "|", pos_bp_str, "|", b1_str, "|", status_str, "|", snps_str, "\n")
}

if (nrow(region_data) > 20) {
  cat("... (", nrow(region_data) - 20, " positions omitted) ...\n")
}

for (i in 1:nrow(last_10)) {
  row <- last_10[i, ]
  pos_10kb_str <- sprintf("%-10.0f", row$pos_10kb)
  pos_bp_str <- sprintf("%-13.0f", row$pos)
  b1_str <- sprintf("%-6.3f", row$B1_freq)
  status_str <- sprintf("%-6s", ifelse(row$estimate_OK == 1, "OK", "FAIL"))
  snps_str <- sprintf("%-4d", row$NSNPs)
  
  cat(pos_10kb_str, "|", pos_bp_str, "|", b1_str, "|", status_str, "|", snps_str, "\n")
}

# Print chromosome statistics
cat("\n=== CHROMOSOME STATISTICS ===\n")
cat("Total positions:", nrow(region_data), "\n")
cat("Positions with estimate_OK = 1:", sum(region_data$estimate_OK == 1, na.rm = TRUE), "\n")
cat("Positions with estimate_OK = 0:", sum(region_data$estimate_OK == 0, na.rm = TRUE), "\n")
cat("Positions with estimate_OK = NA:", sum(is.na(region_data$estimate_OK)), "\n")

# Print SNP count statistics by method
cat("\n=== SNP COUNT STATISTICS ===\n")
snp_stats <- region_data %>%
  group_by(method) %>%
  summarize(
    mean_snps = mean(NSNPs, na.rm = TRUE),
    median_snps = median(NSNPs, na.rm = TRUE),
    min_snps = min(NSNPs, na.rm = TRUE),
    max_snps = max(NSNPs, na.rm = TRUE),
    positions_with_snps = sum(NSNPs > 0, na.rm = TRUE),
    total_positions = n(),
    .groups = "drop"
  )

for (i in 1:nrow(snp_stats)) {
  row <- snp_stats[i, ]
  cat(sprintf("%s: mean=%.1f, median=%.0f, range=[%d-%d], positions with SNPs=%d/%d\n",
              row$method, row$mean_snps, row$median_snps, 
              row$min_snps, row$max_snps, row$positions_with_snps, row$total_positions))
}

cat("\n✓ Chromosome plotting complete!\n")
