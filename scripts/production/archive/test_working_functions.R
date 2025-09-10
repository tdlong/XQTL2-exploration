#!/usr/bin/env Rscript

# Test wrapper for working haplotype estimation functions
# Tests with JUICE data: 1st sample, first 10 positions

library(tidyverse)

# Source the working functions
source("scripts/production/haplotype_estimation_working.R")

# Test parameters
chr <- "chr2R"
method <- "adaptive"
parameter <- 4
output_dir <- "process/JUICE"
param_file <- "helpfiles/JUICE_haplotype_parameters.R"
n_positions <- 100
n_samples <- 1

cat("=== TESTING WORKING HAPLOTYPE ESTIMATION FUNCTIONS ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Testing positions:", n_positions, "\n")
cat("Testing samples:", n_samples, "\n")
cat("================================\n\n")

# 1. Load parameters
cat("1. Loading parameters...\n")
if (!file.exists(param_file)) {
  cat("Error: Parameter file not found:", param_file, "\n")
  quit(status = 1)
}

source(param_file)
cat("✓ Parameter file:", param_file, "\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), samples (", length(names_in_bam), ")\n")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and process data
cat("\n2. Loading and processing data...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  cat("Error: RefAlt file not found:", refalt_file, "\n")
  quit(status = 1)
}

df3 <- process_refalt_data(refalt_file, founders)

# 3. Define test positions
cat("\n3. Defining test positions...\n")
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchrom_start <- euchromatin_boundaries[[chr]][1]
euchrom_end <- euchromatin_boundaries[[chr]][2]

scan_start <- ceiling(euchrom_start / step) * step
scan_end <- floor(euchrom_end / step) * step
all_positions <- seq(scan_start, scan_end, by = step)
test_positions <- head(all_positions, n_positions)

cat("✓ Euchromatin boundaries:", euchrom_start, "to", euchrom_end, "bp\n")
cat("✓ Scan grid:", scan_start, "to", scan_end, "by", step, "\n")
cat("✓ Test positions:", length(test_positions), "\n")
cat("✓ Positions:", paste(test_positions, collapse=", "), "\n")

# 4. Get test samples
test_samples <- head(names_in_bam, n_samples)
cat("✓ Test samples:", paste(test_samples, collapse=", "), "\n")

# 5. Run haplotype estimation
cat("\n4. Running haplotype estimation...\n")
cat("Processing", length(test_positions), "positions ×", length(test_samples), "samples...\n")

results <- expand_grid(
  pos = test_positions,
  sample_name = test_samples
) %>%
  purrr::pmap_dfr(~ {
    cat("\n--- Processing pos:", ..1, "sample:", ..2, "---\n")
    result <- estimate_haplotypes_list_format(
      pos = ..1,
      sample_name = ..2,
      df3 = df3,
      founders = founders,
      h_cutoff = parameter,
      method = method,
      window_size_bp = NULL,
      chr = chr,
      verbose = 1  # Show progress
    )
    
    cat("Result:", ifelse(is.null(result$Haps), "FAILED", "SUCCESS"), "\n")
    if (!is.null(result$Haps)) {
      cat("Haplotype frequencies:", paste(round(result$Haps, 3), collapse=", "), "\n")
      cat("Sum:", round(sum(result$Haps, na.rm = TRUE), 6), "\n")
    }
    
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

# 6. Show results
cat("\n5. Results summary...\n")
cat("✓ Results:", nrow(results), "rows\n")
cat("✓ Positions:", length(unique(results$pos)), "unique positions\n")
cat("✓ Samples:", length(unique(results$sample)), "unique samples\n")

# Show detailed results
cat("\n=== DETAILED RESULTS ===\n")
for (i in 1:nrow(results)) {
  cat("\nPosition:", results$pos[i], "Sample:", results$sample[i], "\n")
  haps <- results$Haps[[i]]
  if (!is.null(haps) && !all(is.na(haps))) {
    cat("  Haplotype frequencies:\n")
    for (j in seq_along(haps)) {
      cat("    ", names(haps)[j], ":", round(haps[j], 4), "\n")
    }
    cat("  Sum:", round(sum(haps, na.rm = TRUE), 6), "\n")
  } else {
    cat("  No haplotype frequencies available\n")
  }
}

cat("\n=== TEST COMPLETE ===\n")
cat("Working functions tested successfully!\n")
