#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/create_summary_file.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== CREATING SUMMARY FILE ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Define methods
methods <- c("fixed_20kb", "fixed_50kb", "fixed_100kb", "fixed_200kb", "fixed_500kb", 
             "adaptive_h4", "adaptive_h6", "adaptive_h8", "adaptive_h10")

# Get all positions divisible by 10kb (10,000)
cat("Loading haplotype data to determine positions...\n")
first_file <- file.path(results_dir, "fixed_20kb_results_chr2R.RDS")
if (!file.exists(first_file)) {
  stop("Haplotype results not found: ", first_file)
}

first_data <- readRDS(first_file)
all_positions <- sort(unique(first_data$pos))
positions_10kb <- all_positions[all_positions %% 10000 == 0]  # Only positions divisible by 10kb

cat("Found", length(positions_10kb), "positions divisible by 10kb\n")
cat("Position range:", min(positions_10kb), "-", max(positions_10kb), "\n\n")

# Initialize summary data frame
summary_data <- data.frame()

# Process each method
for (method in methods) {
  cat("Processing", method, "...\n")
  
  # Extract estimator name
  if (grepl("fixed_", method)) {
    size <- gsub("fixed_", "", method)
    estimator <- paste0("fixed_", size)
    haplo_file <- file.path(results_dir, paste0("fixed_window_", size, "kb_results_", chr, ".RDS"))
  } else {
    h <- gsub("adaptive_h", "", method)
    estimator <- paste0("adaptive_h", h)
    haplo_file <- file.path(results_dir, paste0("adaptive_window_h", h, "_results_", chr, ".RDS"))
  }
  
  snp_file <- file.path(results_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
  
  # Load haplotype data
  if (file.exists(haplo_file)) {
    haplo_data <- readRDS(haplo_file)
    
    # Load SNP data
    if (file.exists(snp_file)) {
      snp_data <- readRDS(snp_file)
      
      # Process each position
      for (pos in positions_10kb) {
        # Get haplotype data for this position
        pos_haplo <- haplo_data %>%
          filter(pos == !!pos) %>%
          select(sample, B1, B2, B3, B4, B5, B6, B7, B8, estimate_OK)
        
        if (nrow(pos_haplo) > 0) {
          # Get SNP data within ±5kb of this position
          pos_snps <- snp_data %>%
            filter(pos >= pos - 5000 & pos <= pos + 5000)
          
          # Calculate RMSE and SNP count for each sample
          for (sample_name in unique(pos_haplo$sample)) {
            sample_haplo <- pos_haplo %>% filter(sample == sample_name)
            sample_snps <- pos_snps %>% filter(sample == sample_name)
            
            if (nrow(sample_snps) > 0) {
              rmse_val <- sqrt(mean((sample_snps$observed - sample_snps$imputed)^2, na.rm = TRUE))
              n_snps <- nrow(sample_snps)
            } else {
              rmse_val <- NA
              n_snps <- 0
            }
            
            # Add to summary data
            summary_row <- data.frame(
              method = method,
              pos = pos,
              pos_10kb = pos / 10000,
              sample = sample_name,
              B1 = sample_haplo$B1,
              B2 = sample_haplo$B2,
              B3 = sample_haplo$B3,
              B4 = sample_haplo$B4,
              B5 = sample_haplo$B5,
              B6 = sample_haplo$B6,
              B7 = sample_haplo$B7,
              B8 = sample_haplo$B8,
              estimate_OK = sample_haplo$estimate_OK,
              n_snps = n_snps,
              rmse = rmse_val
            )
            
            summary_data <- rbind(summary_data, summary_row)
          }
        }
      }
      
      cat("  ✓ Processed", nrow(summary_data %>% filter(method == !!method)), "rows\n")
    } else {
      cat("  ❌ Missing SNP file:", snp_file, "\n")
    }
  } else {
    cat("  ❌ Missing haplotype file:", haplo_file, "\n")
  }
}

# Save summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
saveRDS(summary_data, summary_file)

cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total rows:", nrow(summary_data), "\n")
cat("Unique positions:", length(unique(summary_data$pos)), "\n")
cat("Unique samples:", length(unique(summary_data$sample)), "\n")
cat("Unique methods:", length(unique(summary_data$method)), "\n")

# Show estimate_OK statistics by method
cat("\nEstimate_OK statistics by method:\n")
ok_stats <- summary_data %>%
  group_by(method) %>%
  summarise(
    total = n(),
    ok = sum(estimate_OK == 1, na.rm = TRUE),
    fail = sum(estimate_OK == 0, na.rm = TRUE),
    na_count = sum(is.na(estimate_OK)),
    ok_rate = ok / total,
    .groups = "drop"
  ) %>%
  arrange(desc(ok_rate))

print(ok_stats)

# Show SNP coverage statistics
cat("\nSNP coverage statistics:\n")
snp_stats <- summary_data %>%
  group_by(method) %>%
  summarise(
    total_positions = n(),
    positions_with_snps = sum(n_snps > 0, na.rm = TRUE),
    avg_snps_per_position = mean(n_snps[n_snps > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(positions_with_snps))

print(snp_stats)

cat("\n✓ Summary file saved to:", summary_file, "\n")
cat("File size:", round(file.size(summary_file) / 1024^2, 2), "MB\n")

# Create a quick example of how to use this file
cat("\n=== USAGE EXAMPLE ===\n")
cat("# Load the summary file:\n")
cat("summary_data <- readRDS('", summary_file, "')\n")
cat("\n# Filter to a specific region:\n")
cat("region_data <- summary_data %>% filter(pos >= 5550000 & pos <= 5750000)\n")
cat("\n# Get positions where fixed_20kb fails:\n")
cat("failed_positions <- summary_data %>% filter(method == 'fixed_20kb' & estimate_OK == 0)\n")
