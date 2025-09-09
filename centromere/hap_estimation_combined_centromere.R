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
if (length(args) < 2 || length(args) > 4) {
  cat("Usage: Rscript hap_estimation_combined_centromere.R <param_file> <input_dir> [debug] [debug2]\n")
  cat("Example: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2\n")
  cat("Debug mode: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2 debug\n")
  cat("Debug2 mode: Rscript hap_estimation_combined_centromere.R ../helpfiles/ZINC2_haplotype_parameters.R ../process/ZINC2 debug debug2\n")
  quit(status = 1)
}

param_file <- args[1]
input_dir <- args[2]
debug_mode <- if (length(args) >= 3 && args[3] == "debug") TRUE else FALSE
debug2_mode <- if (length(args) == 4 && args[4] == "debug2") TRUE else FALSE

cat("=== COMBINED CENTROMERE HAPLOTYPE ESTIMATION ===\n")
cat("Input directory:", input_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Debug mode:", if (debug_mode) "ON (first sample only)" else "OFF (all samples)", "\n")
cat("Debug2 mode:", if (debug2_mode) "ON (1500 proximal 2L SNPs)" else "OFF (all 2L SNPs)", "\n")

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

# Process combined chromosomes
for (chr_group in c("2", "3")) {
  cat("\n--- PROCESSING CHROMOSOME GROUP:", chr_group, "---\n")
  
  # Get centromere positions for this chromosome group
  chr_positions <- centromere_positions %>%
    filter(CHROM %in% paste0("chr", chr_group, c("L", "R"))) %>%
    pull(pos)
  
  if (length(chr_positions) == 0) {
    cat("No centromere positions found for chr", chr_group, "- skipping\n")
    next
  }
  
  cat("Combined chr", chr_group, " positions:", length(chr_positions), "\n")
  cat("Position range:", min(chr_positions), "to", max(chr_positions), "\n\n")
  
  # Load and combine RefAlt data from both arms
  combined_data <- tibble()
  for (arm in c("L", "R")) {
    chr <- paste0("chr", chr_group, arm)
    filein <- file.path(input_dir, paste0("RefAlt.", chr, ".txt"))
    if (file.exists(filein)) {
      chr_data <- read_tsv(filein, col_types = cols(.default = "c")) %>%
        mutate(CHROM = chr)
      
      # Apply debug2 mode (1500 proximal SNPs) to chr2L BEFORE combining
      if (debug2_mode && chr == "chr2L") {
        # Get chr2L centromere positions
        chr2L_positions <- centromere_positions %>%
          filter(CHROM == "chr2L") %>%
          pull(pos)
        
        # Limit to 1500 most proximal SNPs (highest positions)
        if (length(chr2L_positions) > 1500) {
          chr2L_positions <- sort(chr2L_positions, decreasing = TRUE)[1:1500]
          cat("DEBUG2 MODE: Limited chr2L to 1500 most proximal SNPs\n")
        }
        
        # Filter chr2L data to these positions
        chr_data <- chr_data %>%
          filter(POS %in% chr2L_positions)
      }
      
      combined_data <- bind_rows(combined_data, chr_data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    cat("No RefAlt data found for chr", chr_group, "\n")
    next
  }
  
  # Renumber positions from 1 to nrows to avoid overlap issues
  combined_data <- combined_data %>%
    mutate(POS = row_number())
  
  # For now, use all the combined data (we'll add centromere filtering later)
  chr_data <- combined_data
  
  cat("Found", nrow(chr_data), "combined positions\n")
  
  if (nrow(chr_data) < 100) {
    cat("WARNING: Very few positions found - may indicate a problem!\n")
  }
  
  # Process the combined data
  cat("Processing combined chr", chr_group, " data...\n")
  
  # TODO: Add the full haplotype estimation logic here
  # This would be similar to the single chromosome processing but with combined data
  
  # Placeholder for now
  chr_results <- tibble(
    CHROM = paste0("chr", chr_group),
    pos = -99,
    sample = "placeholder",
    Groups = list(1:8),
    Haps = list(rep(0.125, 8)),
    Err = list(diag(8)),
    Names = list(founders)
  )
  
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
