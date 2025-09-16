#!/usr/bin/env Rscript

# Debug script to inspect reshaped error matrices and names structure
# Usage: Rscript scripts/ErrMatrix/debug_reshaped_matrices.R <chr> <output_dir> [N]

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  cat("Usage: Rscript scripts/ErrMatrix/debug_reshaped_matrices.R <chr> <output_dir> [N]\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]
N <- if (length(args) == 3) as.numeric(args[3]) else 5

list_dir <- file.path(output_dir, "haplotype_results_list_format")
orig_file <- file.path(list_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
resh_file <- file.path(list_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

cat("Loading files...\n")
orig <- readRDS(orig_file)
resh <- readRDS(resh_file)

# Unnest reshaped
resh_u <- resh %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names))

# Get first N common rows
key_cols <- c("CHROM","pos","sample")
orig_key <- orig %>% select(all_of(key_cols)) %>% mutate(in_orig = TRUE)
resh_key <- resh_u %>% select(all_of(key_cols)) %>% mutate(in_resh = TRUE)
keys <- full_join(orig_key, resh_key, by = key_cols)
common <- keys %>% filter(in_orig == TRUE, in_resh == TRUE) %>% select(all_of(key_cols)) %>%
  arrange(pos, sample) %>% head(N)

cat("Inspecting first", nrow(common), "rows...\n\n")

# Join and inspect
dfc <- common %>%
  left_join(orig, by = key_cols) %>%
  rename(Groups_o = Groups, Haps_o = Haps, Err_o = Err, Names_o = Names) %>%
  left_join(resh_u, by = key_cols) %>%
  rename(Groups_r = Groups, Haps_r = Haps, Err_r = Err, Names_r = Names)

for (i in 1:nrow(dfc)) {
  row <- dfc[i,]
  cat(sprintf("=== Row %d: %s pos=%d sample=%s ===\n", i, row$CHROM, row$pos, row$sample))
  
  # Original
  No <- row$Names_o[[1]]
  Eo <- row$Err_o[[1]]
  cat("Original:\n")
  cat("  Names:", paste(No, collapse=","), "\n")
  cat("  Names class:", class(No), "length:", length(No), "\n")
  cat("  Err class:", class(Eo), "is.matrix:", is.matrix(Eo), "\n")
  if (is.matrix(Eo)) {
    cat("  Err dim:", nrow(Eo), "x", ncol(Eo), "\n")
    cat("  Err rownames:", paste(rownames(Eo), collapse=","), "\n")
    cat("  Err colnames:", paste(colnames(Eo), collapse=","), "\n")
  }
  
  # Reshaped
  Nr <- row$Names_r[[1]]
  Er <- row$Err_r[[1]]
  cat("Reshaped:\n")
  cat("  Names:", paste(Nr, collapse=","), "\n")
  cat("  Names class:", class(Nr), "length:", length(Nr), "\n")
  cat("  Err class:", class(Er), "is.matrix:", is.matrix(Er), "\n")
  if (is.matrix(Er)) {
    cat("  Err dim:", nrow(Er), "x", ncol(Er), "\n")
    cat("  Err rownames:", paste(rownames(Er), collapse=","), "\n")
    cat("  Err colnames:", paste(colnames(Er), collapse=","), "\n")
  }
  
  # Check alignment
  if (is.matrix(Er) && !is.null(rownames(Er)) && !is.null(colnames(Er))) {
    names_in_Er <- all(No %in% rownames(Er)) && all(No %in% colnames(Er))
    cat("  Names alignment check:", names_in_Er, "\n")
    if (!names_in_Er) {
      missing_rows <- setdiff(No, rownames(Er))
      missing_cols <- setdiff(No, colnames(Er))
      if (length(missing_rows) > 0) cat("    Missing from rownames:", paste(missing_rows, collapse=","), "\n")
      if (length(missing_cols) > 0) cat("    Missing from colnames:", paste(missing_cols, collapse=","), "\n")
    }
  } else {
    cat("  Names alignment check: FALSE (not a matrix or missing dimnames)\n")
  }
  
  cat("\n")
}

cat("Done.\n")
