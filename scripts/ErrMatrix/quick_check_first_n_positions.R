#!/usr/bin/env Rscript

# Quick check: run estimator on the first N positions to validate outputs fast
# Usage: Rscript scripts/ErrMatrix/quick_check_first_n_positions.R <chr> <output_dir> <param_file> [N] [sample_name]

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
  library(MASS)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 5) {
  cat("Usage: Rscript scripts/ErrMatrix/quick_check_first_n_positions.R <chr> <output_dir> <param_file> [N] [sample_name]\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]
param_file <- args[3]
N <- if (length(args) >= 4) as.numeric(args[4]) else 10

# Load parameters
source(param_file, local = TRUE)
stopifnot(exists("founders"))

# Choose sample
if (length(args) == 5) {
  sample_name <- args[5]
} else {
  # Default to first sample in names_in_bam if present; else guess from RefAlt later
  sample_name <- if (exists("names_in_bam")) names_in_bam[[1]] else NA_character_
}

# Source BASE_VAR_WIDE.R functions without triggering CLI parsing
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Prepare RefAlt file
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

# Process RefAlt to wide df3 using production logic
cat("Loading and processing RefAlt...\n")
df3 <- process_refalt_data(refalt_file, founders)

# Pick sample if not provided
if (is.na(sample_name)) {
  candidate_samples <- setdiff(names(df3), c("POS", founders))
  if (length(candidate_samples) == 0) stop("No sample columns found in RefAlt data")
  sample_name <- candidate_samples[[1]]
}

# Build all_positions using the same stepping as production
# Use available POS values directly; take first N unique positions
all_positions <- sort(unique(df3$POS))
if (length(all_positions) == 0) stop("No positions found in RefAlt data")
all_positions <- head(all_positions, N)

cat("=== QUICK CHECK ===\n")
cat("Chromosome:", chr, "\n")
cat("Positions to test:", length(all_positions), "\n")
cat("Sample:", sample_name, "\n")
cat("Founders:", paste(founders, collapse=","), "\n\n")

# Use the same max window as estimator (500kb total span)
max_window <- 500000
half_window <- max_window / 2

# h_cutoff comes from adaptive parameter (default 4)
h_cutoff <- 4
method <- "adaptive"

results <- vector("list", length(all_positions))

for (i in seq_along(all_positions)) {
  testing_position <- all_positions[[i]]
  window_start <- max(1, testing_position - half_window)
  window_end <- testing_position + half_window

  df4 <- df3 %>% dplyr::filter(POS >= window_start & POS <= window_end)

  cat(sprintf("[%d/%d] pos=%d, SNPs=%d ", i, length(all_positions), testing_position, nrow(df4)))
  flush.console()

  # Run estimator (minimal verbosity)
  res <- tryCatch({
    est_haps_var(
      testing_position = testing_position,
      sample_name = sample_name,
      df3 = df4,
      founders = founders,
      h_cutoff = h_cutoff,
      method = method,
      chr = chr,
      verbose = 0
    )
  }, error = function(e) e)

  if (inherits(res, "error")) {
    cat("=> ERROR: ", conditionMessage(res), "\n", sep="")
    results[[i]] <- list(pos = testing_position, ok = FALSE, err = conditionMessage(res))
    next
  }

  # Diagnostics
  tr <- tryCatch(sum(diag(res$Err), na.rm = TRUE), error = function(e) NA_real_)
  ng <- tryCatch(length(res$Groups), error = function(e) NA_integer_)
  cat(sprintf("=> groups=%s, trace(Err)=%s\n", ng, format(tr, digits=7)))
  results[[i]] <- list(pos = testing_position, ok = TRUE, groups = ng, trace = tr)
}

# Summary
ok_n <- sum(purrr::map_lgl(results, ~isTRUE(.x$ok)))
cat(sprintf("\nDone. %d/%d positions completed without error.\n", ok_n, length(results)))
