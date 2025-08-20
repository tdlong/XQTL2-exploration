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
    # Filter to sample and add method info
    results <- results %>%
      filter(sample == sample_name) %>%
      mutate(method = paste0("fixed_", size, "kb"))
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
    # Filter to sample and add method info
    results <- results %>%
      filter(sample == sample_name) %>%
      mutate(method = paste0("adaptive_h", h))
    all_results[[length(all_results) + 1]] <- results
    cat("✓ Loaded:", file_name, "\n")
  } else {
    cat("❌ Missing:", file_name, "\n")
  }
}

# Combine all results
combined_results <- bind_rows(all_results)

# Create plot
cat("\nCreating plot...\n")

# Set a colorblind-friendly palette
method_colors <- c(
  "fixed_20kb" = "#E69F00",
  "fixed_50kb" = "#56B4E9", 
  "fixed_100kb" = "#009E73",
  "fixed_200kb" = "#F0E442",
  "fixed_500kb" = "#D55E00",
  "adaptive_h4" = "#0072B2",
  "adaptive_h6" = "#CC79A7",
  "adaptive_h8" = "#000000",
  "adaptive_h10" = "#990099"
)

p <- ggplot(combined_results, aes(x = pos, y = B1, color = method)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE, big.mark = ",")) +
  labs(
    title = paste("B1 Frequencies Across Methods -", sample_name),
    subtitle = paste("Chromosome:", chr),
    x = "Position",
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

# Save full plot
output_file <- file.path(results_dir, paste0("B1_frequencies_", chr, "_", sample_name, ".png"))
ggsave(output_file, p, width = 12, height = 8, dpi = 300)
cat("✓ Full plot saved:", output_file, "\n")

# Create zoomed plot of a region with rapid changes
# Find region with high variance
region_stats <- combined_results %>%
  mutate(region = floor(pos / 100000) * 100000) %>%  # 100kb regions
  group_by(region) %>%
  summarise(var_B1 = var(B1, na.rm = TRUE)) %>%
  arrange(desc(var_B1))

high_var_region <- region_stats$region[1]
zoom_range <- c(high_var_region, high_var_region + 100000)

p_zoom <- ggplot(
  combined_results %>% filter(pos >= zoom_range[1] & pos <= zoom_range[2]),
  aes(x = pos, y = B1, color = method)
) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +  # Add points to see exact positions
  scale_color_manual(values = method_colors) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE, big.mark = ",")) +
  labs(
    title = paste("B1 Frequencies (100kb Zoom) -", sample_name),
    subtitle = paste("Chromosome:", chr, "Region:", format(zoom_range[1], big.mark=","), "-", format(zoom_range[2], big.mark=",")),
    x = "Position",
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
