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

            # Filter data to the region and first sample
            region_data <- summary_data %>%
              filter(sample == first_sample) %>%
              filter(pos >= region_start_bp & pos <= region_end_bp) %>%
              mutate(pos_10kb = pos / 10000)  # Convert to 10kb units for plotting

cat("Data points in region:", nrow(region_data), "\n")

if (nrow(region_data) == 0) {
  stop("No data found in the specified region")
}

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
            p_haplo <- ggplot(region_data, aes(x = pos_10kb, y = B1_freq, color = method)) +
              geom_line(alpha = 0.7) +
              geom_point(size = 1) +
              scale_color_manual(values = method_colors) +
              scale_x_continuous(
                labels = function(x) format(x, scientific = FALSE, big.mark = ","),
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

# Create RMSE plot (bottom panel) - only for positions with SNPs
rmse_data <- region_data %>%
  filter(!is.na(RMSE) & NSNPs > 0)

            if (nrow(rmse_data) > 0) {
              p_rmse <- ggplot(rmse_data, aes(x = pos_10kb, y = RMSE, color = method)) +
                geom_line(alpha = 0.7) +
                geom_point(size = 1) +
                scale_color_manual(values = method_colors) +
                scale_x_continuous(
                  labels = function(x) format(x, scientific = FALSE, big.mark = ","),
                  breaks = seq(round(region_start_bp/10000), round(region_end_bp/10000), by = 2)  # Every 20kb
                ) +
                labs(
                  title = paste("SNP Imputation RMSE -", first_sample),
                  subtitle = paste("Chromosome:", chr, "(averaged over SNPs within ±5kb of each haplotype position)"),
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
  plot_file <- file.path(results_dir, paste0("summary_plot_", chr, "_", first_sample, "_", midpoint_10kb, "_10kb.RDS"))
  ggsave(plot_file, p_combined, width = 12, height = 8, dpi = 300)
  
  cat("✓ Combined plot saved to:", plot_file, "\n")
} else {
  # Save haplotype plot only if no RMSE data
  plot_file <- file.path(results_dir, paste0("summary_plot_", chr, "_", first_sample, "_", midpoint_10kb, "_10kb.RDS"))
  ggsave(plot_file, p_haplo, width = 12, height = 6, dpi = 300)
  
  cat("⚠️ No RMSE data found, saved haplotype plot only\n")
  cat("✓ Plot saved to:", plot_file, "\n")
}

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

cat("\n✓ Plotting complete!\n")
