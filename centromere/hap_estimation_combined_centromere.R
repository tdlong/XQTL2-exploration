#!/usr/bin/env Rscript

# Combined Centromere Haplotype Estimation
# Combines 2L+2R and 3L+3R data for more robust estimation
# Still applies debug2 mode: 1500 most proximal 2L SNPs

# Usage: Rscript hap_estimation_combined_centromere.R <param_file> <input_dir> [debug] [debug2]
# Example: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2
# Debug mode: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2 debug
# Debug2 mode: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2 debug debug2

library(tidyverse)

# Source the centromere-specific function
source("estimate_haplotypes_centromere.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  cat("Usage: Rscript hap_estimation_combined_centromere.R <param_file> <input_dir> [debug]\n")
  cat("Example: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2\n")
  cat("Debug mode: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2 debug\n")
  quit(status = 1)
}

param_file <- args[1]
input_dir <- args[2]
debug_mode <- if (length(args) == 3 && args[3] == "debug") TRUE else FALSE

cat("=== COMBINED CENTROMERE HAPLOTYPE ESTIMATION ===\n")
cat("Input directory:", input_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Debug mode:", if (debug_mode) "ON (first sample only)" else "OFF (all samples)", "\n")

# Load centromere positions
cat("\nLoading centromere positions from sasha_good.rds...\n")
centromere_positions <- readRDS("sasha_good.rds")
cat("✓ Loaded", nrow(centromere_positions), "centromere positions\n")
cat("✓ Chromosomes:", paste(unique(centromere_positions$CHROM), collapse = ", "), "\n\n")

# Load parameters
cat("Loading parameters...\n")
source(param_file)
cat("✓ Parameter file:", param_file, "\n")
cat("✓ Founders:", paste(founders, collapse = ", "), "\n")
cat("✓ Samples:", length(names_in_bam), "\n\n")

# Collect all results in one tibble
all_results <- tibble()

# Process each arm separately first
arm_data <- list()

for (chr in c("chr2L", "chr2R", "chr3L", "chr3R")) {
  cat("\n--- PROCESSING ARM:", chr, "---\n")
  
  # Get centromere positions for this chromosome
  chr_positions <- centromere_positions %>%
    filter(CHROM == chr) %>%
    pull(pos)
  
  if (length(chr_positions) == 0) {
    cat("No centromere positions found for", chr, "- skipping\n")
    next
  }
  
  # Apply 1500 proximal SNPs filter to chr2L (both debug and normal modes)
  if (chr == "chr2L" && length(chr_positions) > 1500) {
    chr_positions <- sort(chr_positions, decreasing = TRUE)[1:1500]
    cat("Limited chr2L to 1500 most proximal SNPs\n")
  }
  
  cat("Sasha positions for", chr, ":", length(chr_positions), "\n")
  cat("Position range:", min(chr_positions), "to", max(chr_positions), "\n")
  
  # Load RefAlt data
  filein <- file.path(input_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(filein)) {
    cat("No RefAlt data found for", chr, "\n")
    next
  }
  
  # Load RefAlt data (copy from working script)
  df <- read.table(filein, header = TRUE) %>%
    mutate(CHROM = chr)
  
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
  
  # Transform data for haplotype estimation (copy from working script)
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
  
  # Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
  founder_wide <- df2 %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  # CRITICAL: Ensure column order matches founders order (pivot_wider can reorder lexically)
  founder_wide <- founder_wide[, c("POS", founders)]
  
  # Verify column order is correct
  cat("Column order after reordering:", paste(colnames(founder_wide)[-1], collapse=", "), "\n")
  cat("Expected order:", paste(founders, collapse=", "), "\n")
  if (!identical(colnames(founder_wide)[-1], founders)) {
    stop("ERROR: Column order mismatch after reordering! This is a critical bug!")
  }
  
  quality_filtered_positions <- founder_wide %>%
    filter(
      if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
    ) %>%
    pull(POS)
  
  # Filter to quality positions and include sample data
  df3 <- df2 %>%
    filter(POS %in% quality_filtered_positions, name %in% c(founders, names_in_bam))
  
  cat("✓ Quality-filtered positions:", length(quality_filtered_positions), "\n")
  
  # Store the processed arm data
  arm_data[[chr]] <- df3
}

# Now combine arms and run haplotype estimation
for (chr_group in c("2", "3")) {
  cat("\n--- COMBINING CHROMOSOME GROUP:", chr_group, "---\n")
  
  # Get the two arms for this group
  arm1 <- paste0("chr", chr_group, "L")
  arm2 <- paste0("chr", chr_group, "R")
  
  if (!(arm1 %in% names(arm_data)) || !(arm2 %in% names(arm_data))) {
    cat("Missing arm data for chr", chr_group, "- skipping\n")
    next
  }
  
  # Rowbind the two arms and renumber positions
  combined_df3 <- bind_rows(arm_data[[arm1]], arm_data[[arm2]]) %>%
    mutate(POS = row_number())  # Renumber positions from 1 to nrows
  
  cat("Combined", arm1, "and", arm2, ":", nrow(combined_df3), "positions\n")
  
  # Run haplotype estimation for each sample
  samples_to_process <- if (debug_mode) names_in_bam[1] else names_in_bam
  
  chr_results <- expand_grid(
    pos = -99,  # Special position for combined
    sample_name = samples_to_process
  ) %>%
    purrr::pmap_dfr(~ {
      result <- estimate_haplotypes_single_window(
        pos = ..1,
        sample_name = ..2,
        df3 = combined_df3,
        founders = founders,
        h_cutoff = h_cutoff,
        window_size_bp = nrow(combined_df3),
        chr = paste0("chr", chr_group)
      )
      
      if (!is.null(result) && !is.null(result$Haps)) {
        cat("✓", ..2, "haplotype frequencies for chr", chr_group, ":\n")
        for (i in seq_along(founders)) {
          cat(sprintf("  %s: %.4f\n", founders[i], result$Haps[i]))
        }
        cat("  Sum:", sprintf("%.6f", sum(result$Haps)), "\n\n")
      } else {
        cat("✗", ..2, "FAILED for chr", chr_group, "\n")
      }
      
      return(tibble(
        CHROM = paste0("chr", chr_group),
        pos = -99,
        sample = ..2,
        Groups = list(result$Groups),
        Haps = list(result$Haps),
        Err = list(result$Err),
        Names = list(result$Names)
      ))
    })
  
  all_results <- bind_rows(all_results, chr_results)
}

# Save all results
if (nrow(all_results) > 0) {
  output_file <- "combined_centromere_all_results.RDS"
  saveRDS(all_results, output_file)
  cat("\n✓ All results saved:", output_file, "\n")
} else {
  cat("\n✗ No successful results\n")
}

cat("\n=== COMBINED CENTROMERE HAPLOTYPE ESTIMATION COMPLETE ===\n")
