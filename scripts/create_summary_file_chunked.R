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

cat("=== CREATING SUMMARY FILE (SIMPLE LOOP APPROACH) ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Define expected files (matching the working plotting script)
fixed_sizes <- c(20, 50, 100, 200, 500)
h_cutoffs <- c(4, 6, 8, 10)

cat("Processing each estimator using working code approach...\n")

all_summaries <- list()

# Process fixed window methods
for (size in fixed_sizes) {
  method <- paste0("fixed_", size, "kb")
  cat("Processing", method, "...\n")
  
  # Haplotype data (exact working code approach)
  h_data <- readRDS(file.path(results_dir, paste0("fixed_window_", size, "kb_results_", chr, ".RDS"))) %>%
    select(chr, pos, estimate_OK, B1, sample) %>%
    mutate(method = method)
  
  # SNP imputation data (exact working code approach)
  snp_file <- file.path(results_dir, paste0("snp_imputation_fixed_", size, "kb_chr", chr, ".RDS"))
  
  if (file.exists(snp_file)) {
    i_data <- readRDS(snp_file) %>%
      mutate(SE = (observed - imputed)^2) %>%
      filter(!is.na(SE)) %>%
      mutate(
        pos_binned = {
          # Calculate breaks (exact working code)
          breaks <- seq(
            from = floor((min(pos, na.rm = TRUE) - 5000) / 10000) * 10000 + 5000,
            to = ceiling((max(pos, na.rm = TRUE) - 5000) / 10000) * 10000 + 15000,
            by = 10000
          )
          # Calculate midpoints for labels
          midpoints <- (breaks[-length(breaks)] + breaks[-1]) / 2
          # Cut with midpoint labels
          cut(pos, breaks = breaks, labels = midpoints, include.lowest = TRUE, right = FALSE)
        }
      ) %>%
      select(chr, pos_binned, SE, sample) %>%
      group_by(chr, pos_binned, sample) %>%
      summarize(RMSE = sqrt(mean(SE)), NSNPs = n(), .groups = "drop") %>%
      rename(pos = pos_binned) %>% 
      mutate(pos = as.numeric(as.character(pos)))
    
    # Join haplotype and SNP data (exact working code approach)
    s_data <- h_data %>% left_join(i_data, by = c("chr", "pos", "sample"))
    
    all_summaries[[length(all_summaries) + 1]] <- s_data
    cat("  ✓ Completed", method, "\n")
  } else {
    cat("  ❌ Missing SNP file:", snp_file, "\n")
  }
}

# Process adaptive window methods
for (h in h_cutoffs) {
  method <- paste0("adaptive_h", h)
  cat("Processing", method, "...\n")
  
  # Haplotype data (exact working code approach)
  h_data <- readRDS(file.path(results_dir, paste0("adaptive_window_h", h, "_results_", chr, ".RDS"))) %>%
    select(chr, pos, estimate_OK, B1, sample) %>%
    mutate(method = method)
  
  # SNP imputation data (exact working code approach)
  snp_file <- file.path(results_dir, paste0("snp_imputation_adaptive_h", h, "_chr", chr, ".RDS"))
  
  if (file.exists(snp_file)) {
    i_data <- readRDS(snp_file) %>%
      mutate(SE = (observed - imputed)^2) %>%
      filter(!is.na(SE)) %>%
      mutate(
        pos_binned = {
          # Calculate breaks (exact working code)
          breaks <- seq(
            from = floor((min(pos, na.rm = TRUE) - 5000) / 10000) * 10000 + 5000,
            to = ceiling((max(pos, na.rm = TRUE) - 5000) / 10000) * 10000 + 15000,
            by = 10000
          )
          # Calculate midpoints for labels
          midpoints <- (breaks[-length(breaks)] + breaks[-1]) / 2
          # Cut with midpoint labels
          cut(pos, breaks = breaks, labels = midpoints, include.lowest = TRUE, right = FALSE)
        }
      ) %>%
      select(chr, pos_binned, SE, sample) %>%
      group_by(chr, pos_binned, sample) %>%
      summarize(RMSE = sqrt(mean(SE)), NSNPs = n(), .groups = "drop") %>%
      rename(pos = pos_binned) %>% 
      mutate(pos = as.numeric(as.character(pos)))
    
    # Join haplotype and SNP data (exact working code approach)
    s_data <- h_data %>% left_join(i_data, by = c("chr", "pos", "sample"))
    
    all_summaries[[length(all_summaries) + 1]] <- s_data
    cat("  ✓ Completed", method, "\n")
  } else {
    cat("  ❌ Missing SNP file:", snp_file, "\n")
  }
}

# Combine all summaries using row bind
final_summary <- bind_rows(all_summaries) %>%
  # Ensure proper column order and names
  select(chr, pos, method, B1_freq = B1, estimate_OK, RMSE, NSNPs, sample)

# Save summary file
summary_file <- file.path(results_dir, paste0("summary_", chr, ".RDS"))
saveRDS(final_summary, summary_file)

cat("\n=== FINAL SUMMARY ===\n")
cat("Total rows:", nrow(final_summary), "\n")
cat("Unique positions:", length(unique(final_summary$pos)), "\n")
cat("Unique methods:", length(unique(final_summary$method)), "\n")
cat("Unique samples:", length(unique(final_summary$sample)), "\n")

cat("\n✓ Summary file saved to:", summary_file, "\n")
cat("File size:", round(file.size(summary_file) / 1024^2, 2), "MB\n")
