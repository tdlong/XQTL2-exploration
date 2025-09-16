#!/usr/bin/env Rscript

# Compare adaptive original vs reshaped adapt_h4 for one chromosome.
# Reports per (pos,sample):
#  - trace(orig), trace(resh), abs diff, Frobenius norm diff of Err
#  - sum of squared differences of Haps (aligned by Names)
#  - flags when Names order differs (scrambling risk)
# Usage: Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R <chr> <output_dir>

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R <chr> <output_dir>\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]

list_dir <- file.path(output_dir, "haplotype_results_list_format")
orig_file <- file.path(list_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
resh_file <- file.path(list_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

stopifnot(file.exists(orig_file), file.exists(resh_file))

cat("Comparing:\n  orig:", orig_file, "\n  resh:", resh_file, "\n\n")

orig <- readRDS(orig_file)
resh <- readRDS(resh_file)

# Reshaped has one row per CHROM/pos with list-of-lists by sample
# Unnest reshaped to rows per sample
resh_u <- resh %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names))

# Keep only present in both
key_cols <- c("CHROM","pos","sample")
orig_key <- orig %>% select(all_of(key_cols)) %>% mutate(in_orig = TRUE)
resh_key <- resh_u %>% select(all_of(key_cols)) %>% mutate(in_resh = TRUE)
keys <- full_join(orig_key, resh_key, by = key_cols)
common <- keys %>% filter(in_orig == TRUE, in_resh == TRUE) %>% select(all_of(key_cols))

cat("Total orig rows:", nrow(orig), " resh rows:", nrow(resh_u), " common:", nrow(common), "\n\n")

# Join payloads for common rows
dfc <- common %>%
  left_join(orig, by = key_cols) %>%
  rename(Groups_o = Groups, Haps_o = Haps, Err_o = Err, Names_o = Names) %>%
  left_join(resh_u, by = key_cols) %>%
  rename(Groups_r = Groups, Haps_r = Haps, Err_r = Err, Names_r = Names)

compute_metrics <- function(row) {
  # Extract
  Eo <- row$Err_o[[1]]
  Er <- row$Err_r[[1]]
  No <- row$Names_o[[1]]
  Nr <- row$Names_r[[1]]
  Ho <- row$Haps_o[[1]]
  Hr <- row$Haps_r[[1]]

  # Defaults
  trace_o <- NA_real_; trace_r <- NA_real_; trace_abs_diff <- NA_real_
  fro_diff <- NA_real_; haps_ssq <- NA_real_; names_match <- NA

  # Names alignment check
  if (is.character(No) && is.character(Nr)) {
    names_match <- identical(No, Nr)
  }

  # Err metrics
  if (is.matrix(Eo) && is.matrix(Er)) {
    # Align Er to No order
    ErA <- Er[No, No, drop = FALSE]
    trace_o <- sum(diag(Eo), na.rm = TRUE)
    trace_r <- sum(diag(ErA), na.rm = TRUE)
    trace_abs_diff <- abs(trace_r - trace_o)
    D <- ErA - Eo
    fro_diff <- sqrt(sum(D^2, na.rm = TRUE))
  }

  # Haps metrics (align by No)
  if (is.numeric(Ho) && is.numeric(Hr)) {
    HrA <- Hr[No]
    haps_ssq <- sum((HrA - Ho)^2, na.rm = TRUE)
  }

  tibble(
    trace_o = trace_o,
    trace_r = trace_r,
    trace_abs_diff = trace_abs_diff,
    fro_err_diff = fro_diff,
    haps_ssq = haps_ssq,
    names_match = names_match
  )
}

metrics <- purrr::pmap_dfr(dfc, compute_metrics)
out <- bind_cols(dfc %>% select(all_of(key_cols)), metrics)

# Summaries
cat("Trace diffs (na removed):\n")
print(
  out %>%
    filter(!is.na(trace_abs_diff)) %>%
    summarise(
      n = n(),
      mean_trace_o = mean(trace_o),
      mean_trace_r = mean(trace_r),
      mean_abs_diff = mean(trace_abs_diff),
      p95_abs_diff = quantile(trace_abs_diff, 0.95),
      max_abs_diff = max(trace_abs_diff)
    )
)

cat("\nFrobenius norm of Err differences (na removed):\n")
print(
  out %>%
    filter(!is.na(fro_err_diff)) %>%
    summarise(
      n = n(),
      mean_fro = mean(fro_err_diff),
      p95_fro = quantile(fro_err_diff, 0.95),
      max_fro = max(fro_err_diff)
    )
)

cat("\nHaps sum of squared differences (na removed):\n")
print(
  out %>%
    filter(!is.na(haps_ssq)) %>%
    summarise(
      n = n(),
      mean_haps_ssq = mean(haps_ssq),
      p95_haps_ssq = quantile(haps_ssq, 0.95),
      max_haps_ssq = max(haps_ssq)
    )
)

cat("\nRows with name order mismatch (possible scrambling):\n")
print(
  out %>% filter(!is.na(names_match) & !names_match) %>% head(20)
)

cat("\nTop 10 by trace_abs_diff:\n")
print(
  out %>% arrange(desc(trace_abs_diff)) %>% head(10)
)

# Write a CSV summary for further inspection
csv_file <- paste0("compare_adapt_vs_reshaped_", chr, ".csv")
readr::write_csv(out, csv_file)
cat("\nSaved per-row comparison to:", csv_file, "\n")


