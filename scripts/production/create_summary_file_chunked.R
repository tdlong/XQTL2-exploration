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

cat("=== CREATING SUMMARY FILE (PROPER SNP COUNTING) ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Define expected files (matching the working plotting script)
fixed_sizes <- c(20, 50, 100, 200, 500)
h_cutoffs <- c(4, 6, 8, 10)

cat("Processing each estimator using proper SNP counting...\n")

all_summaries <- list()

# Process fixed window methods
for (size in fixed_sizes) {
  method <- paste0("fixed_", size, "kb")
  cat("Processing", method, "...\n")
  
  # Haplotype data
  h_data <- readRDS(file.path(results_dir, paste0("fixed_window_", size, "kb_results_", chr, ".RDS"))) %>%
    select(chr, pos, estimate_OK, B1, sample) %>%
    mutate(method = method)
  
  # SNP imputation data
  snp_file <- file.path(results_dir, paste0("snp_imputation_fixed_", size, "kb_", chr, ".RDS"))
  
  if (file.exists(snp_file)) {
    # Load SNP data and count SNPs within the actual window size for each haplotype position
    snp_data <- readRDS(snp_file)
    
    # For fixed windows, count SNPs within ±window_size/2 of each haplotype position
    window_size_bp <- size * 1000
    
    i_data <- h_data %>%
      select(chr, pos, sample) %>%
      left_join(
        snp_data %>%
          group_by(sample) %>%
          group_modify(~ {
            # For each sample, count SNPs within window of each haplotype position
            .x %>%
              mutate(
                # Calculate distance from each haplotype position
                distance_from_pos = abs(pos - .y$pos[1])
              ) %>%
              filter(distance_from_pos <= window_size_bp/2) %>%
              group_by(pos) %>%
              summarize(
                NSNPs = n(),
                RMSE = sqrt(mean((observed - imputed)^2, na.rm = TRUE)),
                .groups = "drop"
              )
          }, .keep = TRUE) %>%
          ungroup(),
        by = c("chr", "pos", "sample")
      ) %>%
      # Fill missing values for positions with no SNPs
      mutate(
        NSNPs = ifelse(is.na(NSNPs), 0, NSNPs),
        RMSE = ifelse(is.na(RMSE), NA, RMSE)
      )
    
    # Join haplotype and SNP data
    s_data <- h_data %>% 
      left_join(i_data %>% select(-chr, -pos), by = "sample") %>%
      mutate(
        NSNPs = ifelse(is.na(NSNPs), 0, NSNPs),
        RMSE = ifelse(is.na(RMSE), NA, RMSE)
      )
    
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
  
  # Haplotype data
  h_data <- readRDS(file.path(results_dir, paste0("adaptive_window_h", h, "_results_", chr, ".RDS"))) %>%
    select(chr, pos, estimate_OK, B1, sample, final_window_size, n_snps) %>%
    mutate(method = method)
  
  # SNP imputation data
  snp_file <- file.path(results_dir, paste0("snp_imputation_adaptive_h", h, "_", chr, ".RDS"))
  
  if (file.exists(snp_file)) {
    # Load SNP data
    snp_data <- readRDS(snp_file)
    
    # For adaptive windows, use the ACTUAL final window size from haplotype results
    # This is the key insight - each position may have used a different window size!
    
    i_data <- h_data %>%
      select(chr, pos, sample, final_window_size) %>%
      left_join(
        snp_data %>%
          group_by(sample) %>%
          group_modify(~ {
            # For each sample, count SNPs within the ACTUAL window size used at each position
            .x %>%
              mutate(
                # Calculate distance from each haplotype position
                distance_from_pos = abs(pos - .y$pos[1])
              ) %>%
              # Use the actual final window size for this specific position
              filter(distance_from_pos <= .y$final_window_size[1]/2) %>%
              group_by(pos) %>%
              summarize(
                NSNPs = n(),
                RMSE = sqrt(mean((observed - imputed)^2, na.rm = TRUE)),
                .groups = "drop"
              )
          }, .keep = TRUE) %>%
          ungroup(),
        by = c("chr", "pos", "sample")
      ) %>%
      # Fill missing values for positions with no SNPs
      mutate(
        NSNPs = ifelse(is.na(NSNPs), 0, NSNPs),
        RMSE = ifelse(is.na(RMSE), NA, RMSE)
      )
    
    # Join haplotype and SNP data
    s_data <- h_data %>% 
      left_join(i_data %>% select(-chr, -pos, -final_window_size), by = "sample") %>%
      mutate(
        NSNPs = ifelse(is.na(NSNPs), 0, NSNPs),
        RMSE = ifelse(is.na(RMSE), NA, RMSE)
      )
    
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
