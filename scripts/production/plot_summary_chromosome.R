#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
              stop("Usage: Rscript scripts/plot_summary_chromosome.R <chr> <param_file> <output_dir>")
}

            chr <- args[1]
            param_file <- args[2]
            output_dir <- args[3]

cat("=== PLOTTING SUMMARY CHROMOSOME ===\n")
            cat("Chromosome:", chr, "\n")
            cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

# Load the summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))

if (!file.exists(summary_file)) {
  stop("Summary file not found: ", summary_file)
}

cat("Loading summary file...\n")
summary_data <- readRDS(summary_file)

            cat("✓ Summary data loaded:", nrow(summary_data), "rows\n")
            cat("Samples available:", paste(unique(summary_data$sample), collapse = ", "), "\n")
            
            # Show position range in data
            pos_range <- range(summary_data$pos)
            cat("Position range in data:", format(pos_range[1], big.mark=","), "to", format(pos_range[2], big.mark=","), "base pairs\n")
            cat("Position range in 10kb units:", round(pos_range[1]/10000), "to", round(pos_range[2]/10000), "\n\n")

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

# Set color scheme for adaptive_h4 only
method_colors <- c("adaptive_h4" = "black")  # Black

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
  geom_point(size = 1.5, na.rm = TRUE, alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.8, span = 0.05) +
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
  geom_point(size = 1.5, na.rm = TRUE, alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.8, span = 0.05) +
                scale_color_manual(values = method_colors) +
                scale_y_continuous(
                  breaks = seq(0.00, 0.15, by = 0.02),
                  labels = sprintf("%.2f", seq(0.00, 0.15, by = 0.02)),
                  limits = c(0.00, 0.15)
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
  geom_point(size = 1.5, na.rm = TRUE, alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.8, span = 0.05) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(pos_range[1]/10000), round(pos_range[2]/10000), by = 100)  # Every 1Mb (100 * 10kb)
  ) +
  # Use log scale for y-axis with reasonable number of labels
  scale_y_log10(
    breaks = c(100, 200, 500, 1000, 2000),
    labels = c("100", "200", "500", "1000", "2000")
  ) +
  labs(
    title = paste("Number of SNPs per Window -", first_sample),
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
plot_file <- file.path(results_dir, paste0("chromosome_summary_", chr, "_", first_sample, "_adaptive_h4.png"))
ggsave(plot_file, p_combined, width = 16, height = 12, dpi = 300)

cat("✓ Combined three-panel plot saved to:", plot_file, "\n")

            # Print summary table for the chromosome
            cat("\n=== SUMMARY TABLE FOR CHROMOSOME ===\n")
            cat("Pos (10kb) | Pos (bp)      | Method        | B1     | Status | SNPs\n")
            cat("------------|---------------|---------------|--------|--------|------\n")

            for (i in 1:nrow(region_data)) {
              row <- region_data[i, ]
              pos_10kb_str <- sprintf("%-10.0f", row$pos_10kb)
              pos_bp_str <- sprintf("%-13.0f", row$pos)
              method_str <- sprintf("%-13s", row$method)
              b1_str <- sprintf("%-6.3f", row$B1_freq)
              status_str <- sprintf("%-6s", ifelse(row$estimate_OK == 1, "OK", "FAIL"))
              snps_str <- sprintf("%-4d", row$NSNPs)

              cat(pos_10kb_str, "|", pos_bp_str, "|", method_str, "|", b1_str, "|", status_str, "|", snps_str, "\n")
            }

# Print region statistics
cat("\n=== REGION STATISTICS ===\n")
cat("Total positions:", nrow(region_data), "\n")
cat("Positions with estimate_OK = 1:", sum(region_data$estimate_OK == 1, na.rm = TRUE), "\n")
cat("Positions with estimate_OK = 0:", sum(region_data$estimate_OK == 0, na.rm = TRUE), "\n")
cat("Positions with estimate_OK = NA:", sum(is.na(region_data$estimate_OK)), "\n")

# Print SNP count statistics by method
cat("\n=== SNP COUNT STATISTICS BY METHOD ===\n")
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

cat("\n✓ Plotting complete!\n")
