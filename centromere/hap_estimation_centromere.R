#!/usr/bin/env Rscriptried 

# CENTROMERE HAPLOTYPE ESTIMATION
# 
# This script runs haplotype estimation for centromere regions
# Processes only chr2L, chr2R, chr3L, chr3R using positions from sasha_good.rds
# One window per chromosome - all SNPs in the subsetted centromere region
#
# Usage: Rscript hap_estimation_centromere.R <param_file> <input_dir>
# Example: Rscript hap_estimation_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2

library(tidyverse)

# Source the centromere-specific function
source("estimate_haplotypes_centromere.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  cat("Usage: Rscript hap_estimation_centromere.R <param_file> <input_dir> [debug]\n")
  cat("Example: Rscript hap_estimation_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2\n")
  cat("Debug mode: Rscript hap_estimation_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2 debug\n")
  quit(status = 1)
}

param_file <- args[1]
input_dir <- args[2]
debug_mode <- if (length(args) == 3 && args[3] == "debug") TRUE else FALSE

# Chromosomes to process
chromosomes <- c("chr2L", "chr2R", "chr3L", "chr3R")

cat("=== CENTROMERE HAPLOTYPE ESTIMATION ===\n")
cat("Input directory:", input_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Debug mode:", if (debug_mode) "ON (first sample only)" else "OFF (all samples)", "\n")
cat("Chromosomes:", paste(chromosomes, collapse = ", "), "\n\n")

# Load centromere positions from sasha_good.rds
cat("Loading centromere positions from sasha_good.rds...\n")
if (!file.exists("sasha_good.rds")) {
  cat("Error: sasha_good.rds not found in current directory\n")
  quit(status = 1)
}

centromere_positions <- readRDS("sasha_good.rds")
cat("✓ Loaded", nrow(centromere_positions), "centromere positions\n")
cat("✓ Chromosomes:", paste(unique(centromere_positions$CHROM), collapse = ", "), "\n\n")

# Load parameters
cat("Loading parameters...\n")
if (!file.exists(param_file)) {
  cat("Error: Parameter file not found:", param_file, "\n")
  quit(status = 1)
}

source(param_file)
cat("✓ Parameter file:", param_file, "\n")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")
cat("✓ Samples:", length(names_in_bam), "\n\n")

# Collect all results in one tibble
all_results <- tibble()

# Process each chromosome
for (chr in chromosomes) {
  cat("\n--- PROCESSING CHROMOSOME:", chr, "---\n")
  
  # Get centromere positions for this chromosome
  chr_positions <- centromere_positions %>%
    filter(CHROM == chr) %>%
    pull(pos)
  
  if (length(chr_positions) == 0) {
    cat("No centromere positions found for", chr, "- skipping\n")
    next
  }
  
  cat("Sasha positions for", chr, ":", length(chr_positions), "\n")
  cat("Position range:", min(chr_positions), "to", max(chr_positions), "\n\n")

  # Load RefAlt data
  cat("Loading RefAlt data...\n")
  filein <- file.path(input_dir, paste0("RefAlt.", chr, ".txt"))
  
  if (!file.exists(filein)) {
    cat("Error: REFALT file not found:", filein, "\n")
    next
  }
  
  df <- read.table(filein, header = TRUE)
  # Check subsetting - how many sasha positions are in RefAlt?
  sasha_in_refalt <- sum(chr_positions %in% df$POS)
  subsetting_fraction <- sasha_in_refalt / length(chr_positions)
  
  cat("✓ Found", sasha_in_refalt, "of", length(chr_positions), "centromere positions (", sprintf("%.1f%%", subsetting_fraction*100), ")\n")
  
  if (subsetting_fraction < 0.8) {
    cat("  WARNING: Low subsetting fraction - may indicate a problem!\n")
  }
  
  # Subset RefAlt to only sasha positions
  df_subset <- df %>%
    filter(POS %in% chr_positions)
  
  # Transform data for haplotype estimation
  df2 <- df_subset %>%
    pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
    mutate(
      RefAlt = str_sub(lab, 1, 3),
      name = str_sub(lab, 5)
    ) %>%
    select(-lab) %>%
    pivot_wider(names_from = RefAlt, values_from = count) %>%
    mutate(
      freq = REF / (REF + ALT),
      N = REF + ALT
    ) %>%
    select(-c("REF", "ALT")) %>%
    as_tibble()
  
  # Filter to founders and samples only
  df3 <- df2 %>%
    filter(name %in% c(founders, names_in_bam))
  
  cat("✓ Using", length(unique(df3$POS)), "positions for haplotype estimation\n")
  
  # Run haplotype estimation - ONE WINDOW per chromosome (all SNPs in subset)
  cat("Running haplotype estimation (one window per chromosome)...\n")
  
  # Use average position rounded to nearest 10kb for the window
  window_center <- round(mean(chr_positions) / 10000) * 10000
  
  # Process samples (debug mode = first sample only, normal mode = all samples)
  samples_to_process <- if (debug_mode) names_in_bam[1] else names_in_bam
  
  if (debug_mode) {
    cat("DEBUG MODE: Processing only first sample:", samples_to_process, "\n")
  } else {
    cat("Processing all", length(samples_to_process), "samples\n")
  }
  
  # Process samples using the same pattern as working production script
  chr_results <- expand_grid(
    pos = window_center,
    sample_name = samples_to_process
  ) %>%
    purrr::pmap_dfr(~ {
      result <- estimate_haplotypes_single_window(
        pos = ..1,
        sample_name = ..2,
        df3 = df3,
        founders = founders,
        h_cutoff = 4,
        window_size_bp = max(chr_positions) - min(chr_positions) + 1,
        chr = chr
      )
      
      # Show useful information: haplotype estimates
      if (!is.null(result) && !is.null(result$Haps)) {
        cat("✓", ..2, "haplotype frequencies:\n")
        for (i in seq_along(founders)) {
          cat(sprintf("  %s: %.4f\n", founders[i], result$Haps[i]))
        }
        cat("  Sum:", sprintf("%.6f", sum(result$Haps)), "\n\n")
      } else {
        cat("✗", ..2, "FAILED\n")
      }
      
      # Create natural row format (same as working production script)
      return(tibble(
        CHROM = chr,
        pos = ..1,
        sample = ..2,
        Groups = list(result$Groups),
        Haps = list(result$Haps),
        Err = list(result$Err),
        Names = list(result$Names)
      ))
    })
  
  # Add chromosome results to main results
  all_results <- bind_rows(all_results, chr_results)
  
  # Clear memory
  rm(df, df_subset, df2, df3)
  gc()
  
}  # End of chromosome loop

# Save all results in one file
if (nrow(all_results) > 0) {
  output_file <- "centromere_all_results.RDS"
  saveRDS(all_results, output_file)
  cat("\n✓ All results saved:", output_file, "\n")
  cat("✓ Total results:", nrow(all_results), "rows\n")
  cat("✓ Chromosomes:", length(unique(all_results$CHROM)), "\n")
  cat("✓ Samples:", length(unique(all_results$sample)), "\n")
} else {
  cat("\n✗ No successful results\n")
}

cat("\n=== CENTROMERE HAPLOTYPE ESTIMATION COMPLETE ===\n")

