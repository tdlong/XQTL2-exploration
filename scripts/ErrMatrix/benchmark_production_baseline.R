#!/usr/bin/env Rscript

# Simple benchmark of production code for first 100 positions
# This gives us a baseline to compare against

library(tidyverse)
library(limSolve)
library(MASS)

# Source the production functions (clean version without main execution)
source("scripts/ErrMatrix/production_functions_only.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript benchmark_production_baseline.R <chr> <output_dir> <param_file>")
}

chr <- args[1]
output_dir <- args[2]
param_file <- args[3]

cat("=== BENCHMARKING PRODUCTION CODE ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Load parameters
source(param_file)

# Load RefAlt data using production function
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
cat("Loading data from:", refalt_file, "\n")

df3 <- process_refalt_data(refalt_file, founders)

# Get sample names
sample_names <- unique(df3$name)
sample_names <- sample_names[!sample_names %in% founders]

cat("Found", length(sample_names), "samples\n")
cat("Found", nrow(df3), "data rows\n\n")

# Define test positions (first 100 positions every 10kb)
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchrom_start <- euchromatin_boundaries[[chr]][1]
euchrom_end <- euchromatin_boundaries[[chr]][2]

# Test positions every 10kb, but only first 100
test_positions <- seq(from = euchrom_start, to = euchrom_end, by = 10000)[1:100]

cat("Testing", length(test_positions), "positions\n")
cat("First position:", test_positions[1], "\n")
cat("Last position:", test_positions[length(test_positions)], "\n\n")

# Time the production function
start_time <- Sys.time()

results <- map_dfr(test_positions, function(test_pos) {
  cat("Processing position:", test_pos, "\n")
  
  # Process first 5 samples only (to keep it fast)
  sample_subset <- sample_names[1:min(5, length(sample_names))]
  
  sample_results <- map_dfr(sample_subset, function(sample_name) {
    tryCatch({
      result <- estimate_haplotypes_list_format(test_pos, sample_name, df3, founders, 4)
      if (!is.null(result)) {
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = sample_name,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      } else {
        return(tibble())
      }
    }, error = function(e) {
      cat("Error processing", sample_name, "at position", test_pos, ":", e$message, "\n")
      return(tibble())
    })
  })
  
  return(sample_results)
})

end_time <- Sys.time()
total_time <- end_time - start_time

cat("\n=== BENCHMARK RESULTS ===\n")
cat("Total time:", round(as.numeric(total_time, units = "secs"), 2), "seconds\n")
cat("Positions processed:", length(test_positions), "\n")
cat("Samples per position:", min(5, length(sample_names)), "\n")
cat("Total function calls:", nrow(results), "\n")
cat("Time per function call:", round(as.numeric(total_time, units = "secs") / nrow(results), 4), "seconds\n")

# Save results
output_file <- file.path(output_dir, paste0("production_baseline_", chr, ".RDS"))
saveRDS(results, output_file)
cat("Results saved to:", output_file, "\n")
