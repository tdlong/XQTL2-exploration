#!/usr/bin/env Rscript

# Plot B1 frequencies across all haplotype estimators to check for suspicious similarity
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 6) {
  stop("Usage: Rscript scripts/plot_B1_frequencies.R <chr> <param_file> <output_dir> <sample_name> [zoom_center_10kb] [zoom_width_10kb]")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
sample_name <- args[4]

# Optional zoom parameters
zoom_center_10kb <- if (length(args) >= 5) as.numeric(args[5]) else NULL
zoom_width_10kb <- if (length(args) >= 6) as.numeric(args[6]) else 20  # Default 200kb window

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

# Check if saved subsetted data exists for this region
subset_file <- file.path(results_dir, paste0("zoomed_data_", chr, "_", sample_name, "_", format(zoom_range_bp[1], big.mark=""), "_", format(zoom_range_bp[2], big.mark=""), ".RDS"))

if (file.exists(subset_file)) {
  cat("✓ Loading saved subsetted data for faster processing...\n")
  saved_data <- readRDS(subset_file)
  zoomed_haplo_data <- saved_data$haplo_data
  zoomed_snp_data <- saved_data$snp_data
  zoom_range_bp <- saved_data$zoom_range
  zoom_range_10kb <- saved_data$zoom_range_10kb
  cat("✓ Loaded data for region:", format(zoom_range_bp[1], big.mark=","), "-", format(zoom_range_bp[2], big.mark=","), "\n")
} else {
  cat("Loading haplotype results...\n")
  
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
  
  # Convert positions to 10kb units for cleaner x-axis
  combined_results <- combined_results %>%
    mutate(pos_10kb = pos / 10000)
  
  # Filter data for zoomed region
  zoomed_haplo_data <- combined_results %>% 
    filter(pos_10kb >= zoom_range_10kb[1] & pos_10kb <= zoom_range_10kb[2])
  
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
    
    # Filter SNP data for zoomed region
    zoomed_snp_data <- snp_combined %>%
      filter(pos_10kb >= zoom_range_10kb[1] & pos_10kb <= zoom_range_10kb[2])
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

# Create zoomed two-panel plot
# Determine zoom region
if (!is.null(zoom_center_10kb)) {
  # Use custom zoom region
  cat("Using custom zoom region: center =", zoom_center_10kb, ", width =", zoom_width_10kb, "\n")
  zoom_start_10kb <- zoom_center_10kb - zoom_width_10kb/2
  zoom_end_10kb <- zoom_center_10kb + zoom_width_10kb/2
  zoom_range_10kb <- c(zoom_start_10kb, zoom_end_10kb)
  zoom_range_bp <- c(zoom_start_10kb * 10000, zoom_end_10kb * 10000)
  zoom_description <- paste("Custom region centered at", zoom_center_10kb)
} else {
  # Find region with high variance
  cat("Finding region with highest variance...\n")
  region_stats <- combined_results %>%
    mutate(region = floor(pos_10kb / 10) * 10) %>%  # 100kb regions (10 units of 10kb)
    group_by(region) %>%
    summarise(var_B1 = var(B1, na.rm = TRUE)) %>%
    arrange(desc(var_B1))

  high_var_region_10kb <- region_stats$region[1]
  zoom_range_10kb <- c(high_var_region_10kb, high_var_region_10kb + zoom_width_10kb)
  zoom_range_bp <- c(high_var_region_10kb * 10000, (high_var_region_10kb + zoom_width_10kb) * 10000)
  zoom_description <- paste("Highest variance region starting at", high_var_region_10kb)
}

# Filter data for zoomed region
zoomed_haplo_data <- combined_results %>% 
  filter(pos_10kb >= zoom_range_10kb[1] & pos_10kb <= zoom_range_10kb[2])

# Create haplotype frequency plot for zoomed region (top panel)
p_zoom_haplo <- ggplot(zoomed_haplo_data, aes(x = pos_10kb, y = B1, color = method)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = ","),
    breaks = seq(zoom_range_10kb[1], zoom_range_10kb[2], by = 2)  # Every 20kb
  ) +
  labs(
    title = paste("B1 Haplotype Frequencies (", zoom_width_10kb*10, "kb Zoom) -", sample_name),
    subtitle = paste("Chromosome:", chr, "Region:", format(zoom_range_bp[1], big.mark=","), "-", format(zoom_range_bp[2], big.mark=","), "\n", zoom_description),
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

# Create RMSE plot for zoomed region (bottom panel)
if (length(snp_results) > 0) {
  # Filter SNP data for zoomed region
  zoomed_snp_data <- snp_combined %>%
    filter(pos_10kb >= zoom_range_10kb[1] & pos_10kb <= zoom_range_10kb[2])
  
  # Save subsetted data for faster debugging
  subset_file <- file.path(results_dir, paste0("zoomed_data_", chr, "_", sample_name, "_", format(zoom_range_bp[1], big.mark=""), "_", format(zoom_range_bp[2], big.mark=""), ".RDS"))
  saveRDS(list(
    haplo_data = zoomed_haplo_data,
    snp_data = zoomed_snp_data,
    zoom_range = zoom_range_bp,
    zoom_range_10kb = zoom_range_10kb
  ), subset_file)
  cat("✓ Saved subsetted data for debugging:", subset_file, "\n")
  
  # Calculate RMSE for each haplotype position (average over SNPs within ±5kb)
  rmse_by_haplo_position <- data.frame()
  
  for (method in unique(zoomed_haplo_data$method)) {
    haplo_positions <- zoomed_haplo_data %>%
      filter(method == !!method) %>%
      pull(pos_10kb)
    
    method_snp_data <- zoomed_snp_data %>%
      filter(method == !!method)
    
    for (haplo_pos in haplo_positions) {
      # Find SNPs within ±5kb of this haplotype position
      haplo_pos_bp <- haplo_pos * 10000
      nearby_snps <- method_snp_data %>%
        filter(pos >= haplo_pos_bp - 5000 & pos <= haplo_pos_bp + 5000)
      
      if (nrow(nearby_snps) > 0) {
        rmse_val <- sqrt(mean((nearby_snps$observed - nearby_snps$imputed)^2, na.rm = TRUE))
        n_snps <- nrow(nearby_snps)
      } else {
        rmse_val <- NA
        n_snps <- 0
      }
      
      rmse_by_haplo_position <- rbind(rmse_by_haplo_position, data.frame(
        method = method,
        pos_10kb = haplo_pos,
        rmse = rmse_val,
        n_snps = n_snps
      ))
    }
  }
  
  cat("✓ Calculated RMSE for", nrow(rmse_by_haplo_position), "haplotype positions\n")
  cat("✓ Average SNPs per position:", round(mean(rmse_by_haplo_position$n_snps, na.rm=TRUE), 1), "\n")
  
  # Create comprehensive table for the zoomed region
  cat("\n=== COMPREHENSIVE TABLE FOR ZOOMED REGION ===\n")
  
  # Get haplotype data for the zoomed region
  zoomed_haplo_table <- zoomed_haplo_data %>%
    select(method, pos_10kb, sample, B1, estimate_OK) %>%
    arrange(pos_10kb, method)
  
  # Get SNP counts for each haplotype position
  snp_counts <- zoomed_snp_data %>%
    group_by(method, pos_10kb) %>%
    summarise(
      n_snps = n(),
      .groups = "drop"
    )
  
  # Combine haplotype and SNP data
  comprehensive_table <- zoomed_haplo_table %>%
    left_join(snp_counts, by = c("method", "pos_10kb")) %>%
    mutate(
      pos_bp = pos_10kb * 10000,
      estimate_status = case_when(
        estimate_OK == 1 ~ "OK",
        estimate_OK == 0 ~ "FAIL",
        is.na(estimate_OK) ~ "NA"
      )
    ) %>%
    select(method, pos_10kb, pos_bp, sample, B1, estimate_status, n_snps) %>%
    arrange(pos_10kb, method)
  
  # Print table
  cat("Position (10kb) | Position (bp)    | Method        | Sample | B1     | Status | SNPs\n")
  cat("----------------|------------------|---------------|--------|--------|--------|------\n")
  
  for (i in 1:nrow(comprehensive_table)) {
    row <- comprehensive_table[i, ]
    pos_10kb_str <- sprintf("%-14.0f", row$pos_10kb)
    pos_bp_str <- sprintf("%-16.0f", row$pos_bp)
    method_str <- sprintf("%-13s", row$method)
    sample_str <- sprintf("%-6s", row$sample)
    b1_str <- sprintf("%-6.3f", row$B1)
    status_str <- sprintf("%-6s", row$estimate_status)
    snps_str <- sprintf("%-4d", row$n_snps)
    
    cat(pos_10kb_str, "|", pos_bp_str, "|", method_str, "|", sample_str, "|", b1_str, "|", status_str, "|", snps_str, "\n")
  }
  
  # Summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  
  # Estimate_OK summary by method
  ok_summary <- comprehensive_table %>%
    group_by(method) %>%
    summarise(
      total_positions = n(),
      ok_count = sum(estimate_status == "OK"),
      fail_count = sum(estimate_status == "FAIL"),
      na_count = sum(estimate_status == "NA"),
      ok_rate = ok_count / total_positions,
      avg_snps = mean(n_snps, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Method         | Positions | OK | FAIL | NA | OK Rate | Avg SNPs\n")
  cat("---------------|-----------|----|------|----|---------|---------\n")
  for (i in 1:nrow(ok_summary)) {
    row <- ok_summary[i, ]
    method_str <- sprintf("%-13s", row$method)
    cat(method_str, "|", sprintf("%-9d", row$total_positions), "|", 
        sprintf("%-2d", row$ok_count), "|", sprintf("%-4d", row$fail_count), "|", 
        sprintf("%-2d", row$na_count), "|", sprintf("%-7.3f", row$ok_rate), "|", 
        sprintf("%-7.1f", row$avg_snps), "\n")
  }
  
  p_zoom_rmse <- ggplot(rmse_by_haplo_position, aes(x = pos_10kb, y = rmse, color = method)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(
      labels = function(x) format(x, scientific = FALSE, big.mark = ","),
      breaks = seq(zoom_range_10kb[1], zoom_range_10kb[2], by = 2)  # Every 20kb
    ) +
    labs(
      title = paste("SNP Imputation RMSE (", zoom_width_10kb*10, "kb Zoom) -", sample_name),
      subtitle = paste("Chromosome:", chr, "Region:", format(zoom_range_bp[1], big.mark=","), "-", format(zoom_range_bp[2], big.mark=","), "\nAveraged over SNPs within ±5kb of each haplotype position"),
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
  
  # Create combined zoomed two-panel plot
  library(patchwork)
  p_zoom_combined <- p_zoom_haplo / p_zoom_rmse + 
    plot_layout(heights = c(1, 1))
  
  # Save combined zoomed plot
  zoom_file <- file.path(results_dir, paste0("B1_frequencies_and_RMSE_", chr, "_", sample_name, "_zoom.png"))
  ggsave(zoom_file, p_zoom_combined, width = 14, height = 10, dpi = 300)
  cat("✓ Combined zoomed plot saved:", zoom_file, "\n")
  
} else {
  # Save haplotype zoom plot only if no SNP data
  zoom_file <- file.path(results_dir, paste0("B1_frequencies_", chr, "_", sample_name, "_zoom.png"))
  ggsave(zoom_file, p_zoom_haplo, width = 12, height = 8, dpi = 300)
  cat("✓ Zoomed haplotype plot saved:", zoom_file, "\n")
}

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
