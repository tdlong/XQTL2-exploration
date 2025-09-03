#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
              stop("Usage: Rscript scripts/plot_summary_region.R <chr> <param_file> <output_dir> <midpoint_10kb>")
}

            chr <- args[1]
            param_file <- args[2]
            output_dir <- args[3]
            midpoint_10kb <- as.numeric(args[4])

cat("=== PLOTTING SUMMARY REGION ===\n")
            cat("Chromosome:", chr, "\n")
            cat("Output directory:", output_dir, "\n")
            cat("Midpoint (10kb):", midpoint_10kb, "\n")
            cat("Region width: 200kb\n\n")

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

            # Filter to the specified region (200kb centered on midpoint)
            region_start_bp <- midpoint_10kb * 10000 - 100000  # 100kb on each side
            region_end_bp <- midpoint_10kb * 10000 + 100000

            cat("Plotting region:", region_start_bp, "-", region_end_bp, "base pairs\n")
            cat("This corresponds to:", round(region_start_bp/10000), "-", round(region_end_bp/10000), "(10kb units)\n\n")
            cat("Note: Input positions are in 10kb units (e.g., 1080 = 10,800,000 base pairs)\n\n")

            # Filter data to the region, first sample, and only the 3 methods of interest
            region_data <- summary_data %>%
              filter(sample == first_sample) %>%
              filter(pos >= region_start_bp & pos <= region_end_bp) %>%
              filter(method %in% c("adaptive_h4", "fixed_20kb", "fixed_100kb")) %>%
              mutate(pos_10kb = pos / 10000)  # Convert to 10kb units for plotting

cat("Data points in region:", nrow(region_data), "\n")

if (nrow(region_data) == 0) {
  stop("No data found in the specified region")
}

# Set color scheme for the 3 methods only
method_colors <- c(
  "adaptive_h4" = "#006400",   # Dark green
  "fixed_20kb" = "#8B0000",    # Dark red
  "fixed_100kb" = "#FF6347"    # Tomato
)

            # Filter data to only include reliable haplotype estimates
            region_data_reliable <- region_data %>%
              filter(estimate_OK == TRUE)

            cat("✓ Filtered to reliable estimates:", nrow(region_data_reliable), "rows\n")
            cat("  Original data:", nrow(region_data), "rows\n")
            cat("  Reliable estimates:", nrow(region_data_reliable), "rows\n")
            cat("  Unreliable estimates:", nrow(region_data) - nrow(region_data_reliable), "rows\n")

            # Create haplotype frequency plot (top panel)
            p_haplo <- ggplot(region_data_reliable, aes(x = pos_10kb, y = B1_freq, color = method)) +
              geom_point(size = 2, na.rm = TRUE) +
              scale_color_manual(values = method_colors) +
              scale_x_continuous(
                labels = function(x) format(x, scientific = FALSE),
                breaks = seq(round(region_start_bp/10000), round(region_end_bp/10000), by = 2)  # Every 20kb
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

# Create MAE plot (middle panel) - only for positions with SNPs
mae_data <- region_data %>%
  filter(!is.na(MAE) & NSNPs > 0)

            if (nrow(mae_data) > 0) {
              p_mae <- ggplot(mae_data, aes(x = pos_10kb, y = MAE, color = method)) +
                geom_point(size = 2, na.rm = TRUE) +
                scale_color_manual(values = method_colors) +
                scale_x_continuous(
                  labels = function(x) format(x, scientific = FALSE),
                  breaks = seq(round(region_start_bp/10000), round(region_end_bp/10000), by = 2)  # Every 20kb
                ) +
                labs(
                  title = paste("SNP Imputation MAE -", first_sample),
                  subtitle = paste("Chromosome:", chr, "(averaged over SNPs within ±5kb of each haplotype position)"),
                  x = NULL,  # No x-axis label for middle panel
                  y = "MAE",
                  color = "Method"
                ) +
                theme_minimal() +
                theme(
                  plot.title = element_text(size = 12, face = "bold"),
                  legend.position = "none",  # Legend only on bottom panel
                  axis.text.x = element_blank()  # No x-axis text for middle panel
                )
            } else {
              # Create empty MAE plot if no data
              p_mae <- ggplot() +
                annotate("text", x = 0.5, y = 0.5, label = "No MAE data available", size = 6) +
                theme_minimal() +
                theme(
                  axis.text = element_blank(),
                  axis.title = element_blank(),
                  panel.grid = element_blank()
                )
            }

# Create SNP count plot (bottom panel) - show all positions including zeros
p_snps <- ggplot(region_data, aes(x = pos_10kb, y = NSNPs, color = method)) +
  geom_point(size = 2, position = position_jitter(width = 1000, seed = 123), na.rm = TRUE) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = seq(round(region_start_bp/10000), round(region_end_bp/10000), by = 2)  # Every 20kb
  ) +
  # Use log scale for y-axis to better show differences between methods
  scale_y_log10(
    breaks = c(1, 5, 10, 25, 50, 100, 250, 500),
    labels = c("1", "5", "10", "25", "50", "100", "250", "500")
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
plot_file <- file.path(results_dir, paste0("summary_plot_", chr, "_", first_sample, "_", midpoint_10kb, "_10kb.png"))
ggsave(plot_file, p_combined, width = 12, height = 10, dpi = 300)

cat("✓ Combined three-panel plot saved to:", plot_file, "\n")

            # Print summary table for the region
            cat("\n=== SUMMARY TABLE FOR REGION ===\n")
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
