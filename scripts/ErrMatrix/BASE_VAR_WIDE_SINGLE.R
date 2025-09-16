#!/usr/bin/env Rscript

# Single-position, no-write debug runner for est_haps_var using production code.
# - Does NOT write any production outputs
# - Runs exactly one position and one sample
# - Prints key diagnostics (trace(Err), groups) to screen only
#
# Usage:
#   Rscript scripts/ErrMatrix/BASE_VAR_WIDE_SINGLE.R <chr> <position> <sample_name> <output_dir> <param_file> [h_cutoff]

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
  library(MASS)
})

print_df4_summary <- function(df4, founders, sample_name) {
  cat("--- df4 summary ---\n")
  cat("rows:", nrow(df4), "\n")
  cat("POS min/max/mean:", min(df4$POS, na.rm=TRUE), max(df4$POS, na.rm=TRUE), round(mean(df4$POS, na.rm=TRUE), 2), "\n")
  cols <- c(founders, sample_name)
  cols <- cols[cols %in% names(df4)]
  stats <- lapply(cols, function(cn) {
    v <- df4[[cn]]
    c(mean=round(mean(v, na.rm=TRUE),6), sd=round(sd(v, na.rm=TRUE),6), na=sum(is.na(v)))
  })
  stats_df <- as.data.frame(do.call(rbind, stats))
  rownames(stats_df) <- cols
  print(stats_df)
  cat("-------------------\n\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 6) {
  cat("Usage: Rscript scripts/ErrMatrix/BASE_VAR_WIDE_SINGLE.R <chr> <position> <sample_name> <output_dir> <param_file> [h_cutoff]\n")
  quit(status = 1)
}

chr <- args[1]
testing_position <- as.numeric(args[2])
sample_name <- args[3]
output_dir <- args[4]
param_file <- args[5]
h_cutoff_cli <- if (length(args) == 6) as.numeric(args[6]) else 4

cat("=== SINGLE-POSITION DEBUG RUN (NO WRITES) ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", testing_position, "\n")
cat("Sample:", sample_name, "\n")
cat("Output dir (for inputs only):", output_dir, "\n")
cat("Param file:", param_file, "\n\n")
cat("Effective h_cutoff (CLI/default):", h_cutoff_cli, "\n\n")

# Source production code (no modification)
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Load parameters (founders, etc.). In adaptive mode we DO NOT read h_cutoff here.
source(param_file, local = TRUE)
stopifnot(exists("founders"))

# Load RefAlt and process df3 (exact production function)
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) stop("RefAlt file not found: ", refalt_file)
df3 <- process_refalt_data(refalt_file, founders)

# Pre-subset df3 to df4 (largest 500kb window)
max_window <- 500000
window_start <- max(1, testing_position - max_window/2)
window_end <- testing_position + max_window/2
df4 <- df3 %>% dplyr::filter(POS >= window_start & POS <= window_end)

cat("Subset stats: rows=", nrow(df4), " POS:", window_start, "-", window_end, "\n", sep="")
cat("Founders:", paste(founders, collapse=","), " h_cutoff:", h_cutoff_cli, "\n\n")

# Print simple summary (inputs)
print_df4_summary(df4, founders, sample_name)

# Run estimator once
res <- est_haps_var(
  testing_position = testing_position,
  sample_name = sample_name,
  df3 = df4,
  founders = founders,
  h_cutoff = h_cutoff_cli,
  method = "adaptive",
  chr = chr,
  verbose = 1
)

trace_Err <- tryCatch(sum(diag(res$Err), na.rm=TRUE), error=function(e) NA)
cat("=== RESULT ===\n")
cat("Groups:", paste(res$Groups, collapse=","), "\n")
cat("Names:", paste(res$Names, collapse=","), "\n")
cat("trace(Err):", format(trace_Err, digits=8), "\n")
cat("Done.\n")


