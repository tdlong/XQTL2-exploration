#!/usr/bin/env Rscript

# Compare adaptive original vs reshaped adapt_h4 for one chromosome.
# Reports per (pos,sample):
#  - trace(orig), trace(resh), abs diff, Frobenius norm diff of Err
#  - sum of squared differences of Haps (aligned by Names)
#  - flags when Names order differs (scrambling risk)
# Usage:
#   Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R <chr> <output_dir> [limit]

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  cat("Usage: Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R <chr> <output_dir> [limit]\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]
limit_n <- if (length(args) == 3) as.integer(args[3]) else NA_integer_

list_dir <- file.path(output_dir, "haplotype_results_list_format")
orig_file <- file.path(list_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
resh_file <- file.path(list_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

stopifnot(file.exists(orig_file), file.exists(resh_file))

cat("Comparing:\n  orig:", orig_file, "\n  resh:", resh_file, "\n\n")

orig <- readRDS(orig_file)
resh <- readRDS(resh_file)

# Reshaped has one row per CHROM/pos with list-of-lists by sample
# Unnest reshaped to rows per sample
dash_u <- resh %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names))

# Keep only present in both
key_cols <- c("CHROM","pos","sample")
orig_key <- orig %>% select(all_of(key_cols)) %>% mutate(in_orig = TRUE)
dash_key <- dash_u %>% select(all_of(key_cols)) %>% mutate(in_resh = TRUE)
keys <- full_join(orig_key, dash_key, by = key_cols)
common <- keys %>% filter(in_orig == TRUE, in_resh == TRUE) %>% select(all_of(key_cols))

# Apply optional limit for speed
common <- common %>% arrange(pos, sample)
if (!is.na(limit_n) && limit_n > 0) {
  common <- head(common, limit_n)
}

cat("Total orig rows:", nrow(orig), " resh rows:", nrow(dash_u), " common (after limit):", nrow(common), "\n\n")

# Join payloads for common rows
dfc <- common %>%
  left_join(orig, by = key_cols) %>%
  rename(Groups_o = Groups, Haps_o = Haps, Err_o = Err, Names_o = Names) %>%
  left_join(dash_u, by = key_cols) %>%
  rename(Groups_r = Groups, Haps_r = Haps, Err_r = Err, Names_r = Names)

compute_metrics <- function(row) {
  tryCatch({
    # Extract - handle nested list structure correctly
    Eo <- row$Err_o[[1]]
    Er <- row$Err_r[[1]]
    No <- row$Names_o[[1]]
    Nr <- row$Names_r[[1]]
    Ho <- row$Haps_o[[1]]
    Hr <- row$Haps_r[[1]]
    
    # Debug: check if we're getting the right structure
    if (row$pos == 4560000 && row$sample == "Rep01_W_F") {
      cat("DEBUG: No class:", class(No), "length:", length(No), "content:", paste(No, collapse=","), "\n")
      cat("DEBUG: Er rownames length:", length(rownames(Er)), "colnames length:", length(colnames(Er)), "\n")
    }

    # Defaults
    trace_o <- NA_real_; trace_r <- NA_real_; trace_abs_diff <- NA_real_
    fro_diff <- NA_real_; haps_ssq <- NA_real_; names_match <- NA
    has_err_mats <- FALSE; has_names <- FALSE; names_in_Er <- FALSE

    # Names alignment check
    if (!is.null(No)) No <- as.character(No)
    if (!is.null(Nr)) Nr <- as.character(Nr)
    if (is.character(No) && is.character(Nr)) {
      names_match <- identical(No, Nr)
      has_names <- TRUE
    }

    # Err metrics (coerce to matrices and verify names exist)
    if (!is.null(Eo)) Eo <- tryCatch(as.matrix(Eo), error=function(e) NULL)
    if (!is.null(Er)) Er <- tryCatch(as.matrix(Er), error=function(e) NULL)
    if (is.matrix(Eo) && is.matrix(Er) && length(No) > 0) {
      rn <- rownames(Er); cn <- colnames(Er)
      names_in_Er <- !is.null(rn) && !is.null(cn) && all(No %in% rn) && all(No %in% cn)
      # Debug: print first row details
      if (row$pos == 4560000 && row$sample == "Rep01_W_F") {
        cat("DEBUG: Eo is.matrix:", is.matrix(Eo), "Er is.matrix:", is.matrix(Er), "\n")
        cat("DEBUG: No length:", length(No), "rn length:", length(rn), "cn length:", length(cn), "\n")
        cat("DEBUG: names_in_Er:", names_in_Er, "\n")
      }
      if (names_in_Er) {
        ErA <- Er[No, No, drop = FALSE]
        trace_o <- sum(diag(Eo), na.rm = TRUE)
        trace_r <- sum(diag(ErA), na.rm = TRUE)
        trace_abs_diff <- abs(trace_r - trace_o)
        D <- ErA - Eo
        fro_diff <- sqrt(sum(D^2, na.rm = TRUE))
        has_err_mats <- TRUE
      }
    }

    # Haps metrics (align by No when partially available)
    if (!is.null(Ho) && !is.null(Hr)) {
      if (is.list(Hr)) Hr <- unlist(Hr)
      if (is.list(Ho)) Ho <- unlist(Ho)
      if (is.numeric(Ho) && is.numeric(Hr) && length(No) > 0 && !is.null(names(Hr))) {
        names(Hr) <- as.character(names(Hr))
        ok_idx <- No[No %in% names(Hr)]
        if (length(ok_idx) > 0) {
          HrA <- Hr[ok_idx]
          HoA <- Ho[ok_idx]
          haps_ssq <- sum((HrA - HoA)^2, na.rm = TRUE)
        }
      }
    }

    tibble::tibble(
      trace_o = trace_o,
      trace_r = trace_r,
      trace_abs_diff = trace_abs_diff,
      fro_err_diff = fro_diff,
      haps_ssq = haps_ssq,
      names_match = names_match,
      has_err_mats = has_err_mats,
      has_names = has_names,
      names_in_Er = names_in_Er
    )
  }, error = function(e) {
    tibble::tibble(
      trace_o = NA_real_,
      trace_r = NA_real_,
      trace_abs_diff = NA_real_,
      fro_err_diff = NA_real_,
      haps_ssq = NA_real_,
      names_match = NA,
      has_err_mats = FALSE,
      has_names = FALSE,
      names_in_Er = FALSE
    )
  })
}

# Defensive rowwise mapping to ensure single bad row doesn't abort job
metrics <- dfc %>% dplyr::rowwise() %>% dplyr::do(compute_metrics(.)) %>% dplyr::ungroup()
out <- bind_cols(dfc %>% select(all_of(key_cols)), metrics)

valid_rows <- sum(!is.na(out$trace_abs_diff) | !is.na(out$fro_err_diff) | !is.na(out$haps_ssq))
cat("Valid metric rows:", valid_rows, "of", nrow(out), "\n")
if (valid_rows == 0) {
  cat("No valid metrics computed. Diagnostics (first 5 rows):\n")
  print(out %>% select(CHROM, pos, sample, has_err_mats, has_names, names_in_Er) %>% head(5))
}

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

# Write a CSV summary for further inspection (includes diagnostics columns)
csv_file <- paste0("compare_adapt_vs_reshaped_", chr, ".csv")
readr::write_csv(out, csv_file)
cat("\nSaved per-row comparison to:", csv_file, "\n")


