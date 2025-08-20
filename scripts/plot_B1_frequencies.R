#!/usr/bin/env Rscript

# Plot B1 frequencies across all haplotype estimators to check for suspicious similarity
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript plot_B1_frequencies.R <chr> <param_file> <output_dir> <sample_name>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
sample_name <- args[4]

cat("=== PLOTTING B1 FREQUENCIES ACROSS METHODS ===\n")
cat("Chromosome:", chr, "\n")
cat("Sample:", sample_name, "\n")

# Load parameters
source(param_file)
cat("✓ Parameter file loaded\n")

# Define expected files
fixed_sizes <- c(20, 50, 100, 200, 500)
h_cutoffs <- c(4, 6, 8, 10)

results_dir <- file.path(output_dir, "haplotype_results")

# Load all haplotype files
all_results <- list()

# Load fixed window results
for (size in fixed_sizes) {
  file_name <- paste0("fixed_window_", size, "kb_results_", chr, ".RDS")
  file_path <- file.path(results_dir, file_name)
  
  if (file.exists(file_path)) {
    results <- readRDS(file_path)
    
    # Debug: Show raw data structure
    cat("\nDEBUG - Raw data for", file_name, ":\n")
    cat("File size:", file.size(file_path), "bytes\n")
    cat("Columns:", paste(names(results), collapse=", "), "\n")
    cat("First few rows:\n")
    print(head(results))
    
    # Check if this file has unique content
    file_hash <- digest::digest(results, algo="md5")
    cat("File content hash:", substr(file_hash, 1, 8), "\n")
    
    # Filter to sample and add method info
    results <- results %>%
      filter(sample == sample_name) %>%
      mutate(method = paste0("fixed_", size, "kb"))
    
    # Debug: Show filtered data
    cat("\nAfter filtering to", sample_name, ":\n")
    cat("Rows:", nrow(results), "\n")
    cat("Unique positions:", length(unique(results$pos)), "\n")
    cat("Position range:", min(results$pos), "-", max(results$pos), "\n")
    cat("B1 range:", min(results$B1, na.rm=TRUE), "-", max(results$B1, na.rm=TRUE), "\n")
    
    all_results[[length(all_results) + 1]] <- results
    cat("✓ Loaded:", file_name, "\n")
  } else {
    cat("❌ Missing:", file_name, "\n")
  }
}

# Load adaptive window results
for (h in h_cutoffs) {
  file_name <- paste0("adaptive_window_h", h, "_results_", chr, ".RDS")
  file_path <- file.path(results_dir, file_name)
  
  if (file.exists(file_path)) {
    results <- readRDS(file_path)
    
    # Debug: Show raw data structure
    cat("\nDEBUG - Raw data for", file_name, ":\n")
    cat("File size:", file.size(file_path), "bytes\n")
    cat("Columns:", paste(names(results), collapse=", "), "\n")
    cat("First few rows:\n")
    print(head(results))
    
    # Check if this file has unique content
    file_hash <- digest::digest(results, algo="md5")
    cat("File content hash:", substr(file_hash, 1, 8), "\n")
    
    # Filter to sample and add method info
    results <- results %>%
      filter(sample == sample_name) %>%
      mutate(method = paste0("adaptive_h", h))
    
    # Debug: Show filtered data
    cat("\nAfter filtering to", sample_name, ":\n")
    cat("Rows:", nrow(results), "\n")
    cat("Unique positions:", length(unique(results$pos)), "\n")
    cat("Position range:", min(results$pos), "-", max(results$pos), "\n")
    cat("B1 range:", min(results$B1, na.rm=TRUE), "-", max(results$B1, na.rm=TRUE), "\n")
    
    all_results[[length(all_results) + 1]] <- results
    cat("✓ Loaded:", file_name, "\n")
  } else {
    cat("❌ Missing:", file_name, "\n")
  }
}

# Combine all results
combined_results <- bind_rows(all_results)

# Check for identical positions across methods
cat("\n=== CHECKING FOR IDENTICAL POSITIONS ===\n")
positions_by_method <- combined_results %>%
  group_by(method) %>%
  summarise(
    positions = list(sort(unique(pos))),
    n_pos = length(unique(pos))
  )

# Compare each method's positions to fixed_500kb
base_positions <- positions_by_method$positions[positions_by_method$method == "fixed_500kb"][[1]]
for (i in 1:nrow(positions_by_method)) {
  method <- positions_by_method$method[i]
  if (method != "fixed_500kb") {
    current_positions <- positions_by_method$positions[[i]]
    identical_pos <- identical(current_positions, base_positions)
    cat(sprintf("%s: %d positions, identical to fixed_500kb: %s\n", 
                method, positions_by_method$n_pos[i], 
                ifelse(identical_pos, "YES ⚠️", "no")))
  }
}

# Create plot
cat("\nCreating plot...\n")

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

# Convert positions to 10kb units for cleaner x-axis
combined_results <- combined_results %>%
  mutate(pos_10kb = pos / 10000)

# Create haplotype frequency plot (top panel)
p_haplo <- ggplot(combined_results, aes(x = pos_10kb, y = B1, color = method)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = ","),
    breaks = seq(0, max(combined_results$pos_10kb), by = 50)  # Every 500kb
  ) +
  labs(
    title = paste("B1 Haplotype Frequencies -", sample_name),
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

# Load SNP imputation results for bottom panel
cat("\nLoading SNP imputation results for RMSE analysis...\n")
snp_results <- list()

for (method in unique(combined_results$method)) {
  # Extract estimator name from method
  if (grepl("fixed_", method)) {
    size <- gsub("fixed_", "", method)
    estimator <- paste0("fixed_", size)
  } else {
    h <- gsub("adaptive_h", "", method)
    estimator <- paste0("adaptive_h", h)
  }
  
  snp_file <- file.path(results_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
  
  if (file.exists(snp_file)) {
    snp_data <- readRDS(snp_file) %>%
      filter(sample == sample_name) %>%
      mutate(method = method)
    snp_results[[length(snp_results) + 1]] <- snp_data
    cat("✓ Loaded SNP data for", method, "\n")
  } else {
    cat("❌ Missing SNP file for", method, ":", snp_file, "\n")
  }
}

# Combine SNP results and calculate RMSE by position
if (length(snp_results) > 0) {
  snp_combined <- bind_rows(snp_results) %>%
    mutate(pos_10kb = pos / 10000)
  
  # Calculate RMSE for each haplotype position (average over SNPs within ±5kb)
  rmse_by_position <- snp_combined %>%
    group_by(method, pos_10kb) %>%
    summarise(
      rmse = sqrt(mean((observed - imputed)^2, na.rm = TRUE)),
      n_snps = n(),
      .groups = "drop"
    )
  
  # Create RMSE plot (bottom panel)
  p_rmse <- ggplot(rmse_by_position, aes(x = pos_10kb, y = rmse, color = method)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(
      labels = function(x) format(x, scientific = FALSE, big.mark = ","),
      breaks = seq(0, max(rmse_by_position$pos_10kb), by = 50)  # Every 500kb
    ) +
    labs(
      title = paste("SNP Imputation RMSE -", sample_name),
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
  library(patchwork)
  p_combined <- p_haplo / p_rmse + 
    plot_layout(heights = c(1, 1))
  
  # Save combined plot
  output_file <- file.path(results_dir, paste0("B1_frequencies_and_RMSE_", chr, "_", sample_name, ".png"))
  ggsave(output_file, p_combined, width = 14, height = 10, dpi = 300)
  cat("✓ Combined plot saved:", output_file, "\n")
  
} else {
  cat("❌ No SNP imputation data found - saving haplotype plot only\n")
  # Save haplotype plot only
  output_file <- file.path(results_dir, paste0("B1_frequencies_", chr, "_", sample_name, ".png"))
  ggsave(output_file, p_haplo, width = 12, height = 8, dpi = 300)
  cat("✓ Haplotype plot saved:", output_file, "\n")
}

# Create zoomed plot of a region with rapid changes (2x larger window)
# Find region with high variance
region_stats <- combined_results %>%
  mutate(region = floor(pos_10kb / 10) * 10) %>%  # 100kb regions (10 units of 10kb)
  group_by(region) %>%
  summarise(var_B1 = var(B1, na.rm = TRUE)) %>%
  arrange(desc(var_B1))

high_var_region_10kb <- region_stats$region[1]
zoom_range_10kb <- c(high_var_region_10kb, high_var_region_10kb + 20)  # 200kb window (20 units of 10kb)
zoom_range_bp <- c(high_var_region_10kb * 10000, (high_var_region_10kb + 20) * 10000)

p_zoom <- ggplot(
  combined_results %>% filter(pos_10kb >= zoom_range_10kb[1] & pos_10kb <= zoom_range_10kb[2]),
  aes(x = pos_10kb, y = B1, color = method)
) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +  # Add points to see exact positions
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = ","),
    breaks = seq(zoom_range_10kb[1], zoom_range_10kb[2], by = 2)  # Every 20kb
  ) +
  labs(
    title = paste("B1 Frequencies (200kb Zoom) -", sample_name),
    subtitle = paste("Chromosome:", chr, "Region:", format(zoom_range_bp[1], big.mark=","), "-", format(zoom_range_bp[2], big.mark=",")),
    x = "Position (10kb units)",
    y = "B1 Frequency",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save zoomed plot
zoom_file <- file.path(results_dir, paste0("B1_frequencies_", chr, "_", sample_name, "_zoom.png"))
ggsave(zoom_file, p_zoom, width = 12, height = 8, dpi = 300)
cat("✓ Zoomed plot saved:", zoom_file, "\n")

# Print some summary statistics
cat("\nSummary Statistics:\n")

# Verify sample filtering worked
sample_counts <- combined_results %>%
  group_by(method) %>%
  summarise(
    unique_samples = n_distinct(sample),
    sample_list = paste(unique(sample), collapse=", ")
  )
cat("\nSample check:\n")
print(sample_counts)

# Check position spacing
pos_spacing <- combined_results %>%
  group_by(method) %>%
  arrange(pos) %>%
  mutate(pos_diff = pos - lag(pos)) %>%
  summarise(
    min_spacing = min(pos_diff, na.rm = TRUE),
    median_spacing = median(pos_diff, na.rm = TRUE),
    max_spacing = max(pos_diff, na.rm = TRUE),
    total_positions = n()
  )
cat("\nPosition spacing (bp):\n")
print(pos_spacing)

# Basic B1 frequency stats
summary_stats <- combined_results %>%
  group_by(method) %>%
  summarise(
    mean_B1 = mean(B1, na.rm = TRUE),
    sd_B1 = sd(B1, na.rm = TRUE),
    min_B1 = min(B1, na.rm = TRUE),
    max_B1 = max(B1, na.rm = TRUE),
    na_count = sum(is.na(B1)),
    total_pos = n(),
    rapid_changes = sum(abs(B1 - lag(B1)) > 0.5, na.rm = TRUE)  # Count big jumps
  ) %>%
  mutate(
    na_pct = na_count / total_pos * 100,
    rapid_change_pct = rapid_changes / total_pos * 100
  )

cat("\nB1 frequency statistics:\n")
print(summary_stats)

# Check correlation between methods
cat("\nCorrelation between fixed_500kb and other methods:\n")
base_method <- combined_results %>%
  filter(method == "fixed_500kb") %>%
  select(pos, B1) %>%
  rename(base_B1 = B1)

for (m in unique(combined_results$method)) {
  if (m != "fixed_500kb") {
    method_data <- combined_results %>%
      filter(method == m) %>%
      select(pos, B1) %>%
      inner_join(base_method, by = "pos")
    
    cor_val <- cor(method_data$B1, method_data$base_B1, use = "complete.obs")
    cat(sprintf("%s: r = %.3f\n", m, cor_val))
  }
}

# Identify which method is "smooth" vs "oscillating"
cat("\n=== SMOOTHNESS ANALYSIS ===\n")
smoothness_analysis <- combined_results %>%
  group_by(method) %>%
  arrange(pos) %>%
  mutate(
    b1_diff = abs(B1 - lag(B1)),
    b1_diff_squared = b1_diff^2
  ) %>%
  summarise(
    mean_step_change = mean(b1_diff, na.rm = TRUE),
    total_variation = sum(b1_diff_squared, na.rm = TRUE),
    max_step_change = max(b1_diff, na.rm = TRUE)
  ) %>%
  arrange(mean_step_change)

cat("Methods ordered by smoothness (lowest mean step change = smoothest):\n")
for (i in 1:nrow(smoothness_analysis)) {
  method <- smoothness_analysis$method[i]
  mean_change <- smoothness_analysis$mean_step_change[i]
  max_change <- smoothness_analysis$max_step_change[i]
  cat(sprintf("%d. %s: mean=%.4f, max=%.4f\n", i, method, mean_change, max_change))
}
