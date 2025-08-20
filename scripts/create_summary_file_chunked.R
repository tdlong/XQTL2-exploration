#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/create_summary_file_chunked.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== CREATING SIMPLE SUMMARY FILE ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Define expected files (matching the working plotting script)
fixed_sizes <- c(20, 50, 100, 200, 500)
h_cutoffs <- c(4, 6, 8, 10)

cat("Loading haplotype results...\n")

# Load all haplotype files (matching the working plotting script logic)
all_results <- list()

# Load fixed window results
for (size in fixed_sizes) {
  file_name <- paste0("fixed_window_", size, "kb_results_", chr, ".RDS")
  file_path <- file.path(results_dir, file_name)
  
  if (file.exists(file_path)) {
    results <- readRDS(file_path)
    
    # Add method info (matching the working plotting script)
    results <- results %>%
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
    
    # Add method info (matching the working plotting script)
    results <- results %>%
      mutate(method = paste0("adaptive_h", h))
    
    all_results[[length(all_results) + 1]] <- results
    cat("✓ Loaded:", file_name, "\n")
  } else {
    cat("❌ Missing:", file_name, "\n")
  }
}

# Combine all haplotype results
combined_results <- bind_rows(all_results)

# Get positions every 10kb from haplotype files
positions_10kb <- sort(unique(combined_results$pos))
positions_10kb <- positions_10kb[positions_10kb %% 10000 == 0]  # Only positions divisible by 10kb

cat("Found", length(positions_10kb), "positions divisible by 10kb\n")
cat("Position range:", min(positions_10kb), "-", max(positions_10kb), "\n\n")

# Initialize summary data
summary_data <- data.frame()

# Process each position
cat("Processing positions...\n")
for (pos in positions_10kb) {
  cat("  Position:", format(pos, big.mark=","), "\n")
  
  # Get haplotype data for this position
  pos_haplo <- combined_results %>%
    filter(pos == !!pos)
  
  if (nrow(pos_haplo) > 0) {
    # Process each method for this position
    for (method in unique(pos_haplo$method)) {
      # Get haplotype data for this method and position
      method_haplo <- pos_haplo %>%
        filter(method == !!method)
      
      if (nrow(method_haplo) > 0) {
        # Extract estimator name from method (matching the working plotting script)
        if (grepl("fixed_", method)) {
          size <- gsub("fixed_", "", method)
          estimator <- paste0("fixed_", size)
        } else {
          h <- gsub("adaptive_h", "", method)
          estimator <- paste0("adaptive_h", h)
        }
        
        # Load SNP data
        snp_file <- file.path(results_dir, paste0("snp_imputation_", estimator, "_", chr, ".RDS"))
        
        if (file.exists(snp_file)) {
          snp_data <- readRDS(snp_file)
          
          # Get SNPs within ±5kb of this position
          nearby_snps <- snp_data %>%
            filter(pos >= pos - 5000 & pos <= pos + 5000)
          
          if (nrow(nearby_snps) > 0) {
            # Calculate RMSE across all samples
            rmse_val <- sqrt(mean((nearby_snps$observed - nearby_snps$imputed)^2, na.rm = TRUE))
            n_snps <- nrow(nearby_snps)
          } else {
            rmse_val <- NA
            n_snps <- 0
          }
          
          # Add row for each sample
          for (sample_name in unique(method_haplo$sample)) {
            sample_haplo <- method_haplo %>%
              filter(sample == sample_name)
            
            summary_row <- data.frame(
              chr = chr,
              pos = pos,
              method = method,
              B1_freq = sample_haplo$B1,
              estimate_OK = sample_haplo$estimate_OK,
              RMSE = rmse_val,
              NSNPs = n_snps
            )
            
            summary_data <- rbind(summary_data, summary_row)
          }
        }
      }
    }
  }
}

# Save summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
saveRDS(summary_data, summary_file)

cat("\n=== FINAL SUMMARY ===\n")
cat("Total rows:", nrow(summary_data), "\n")
cat("Unique positions:", length(unique(summary_data$pos)), "\n")
cat("Unique methods:", length(unique(summary_data$method)), "\n")
cat("Unique samples:", length(unique(summary_data$sample)), "\n")

cat("\n✓ Summary file saved to:", summary_file, "\n")
cat("File size:", round(file.size(summary_file) / 1024^2, 2), "MB\n")
