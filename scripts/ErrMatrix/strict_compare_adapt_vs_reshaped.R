#!/usr/bin/env Rscript

# Strict per (CHROM,pos,sample) equality check between adaptive original and reshaped adapt_h4.
# Pass/fail accounting with tight tolerance. Writes only mismatches CSV.
#
# Usage:
#   Rscript scripts/ErrMatrix/strict_compare_adapt_vs_reshaped.R <chr> <output_dir> [tolerance]

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  cat("Usage: Rscript scripts/ErrMatrix/strict_compare_adapt_vs_reshaped.R <chr> <output_dir> [tolerance]\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]
tol <- if (length(args) == 3) as.numeric(args[3]) else 1e-10

base_dir <- file.path(output_dir, "haplotype_results_list_format")
orig_file <- file.path(base_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
resh_file <- file.path(base_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

stopifnot(file.exists(orig_file), file.exists(resh_file))

orig <- readRDS(orig_file)
resh <- readRDS(resh_file) %>% tidyr::unnest(c(sample, Groups, Haps, Err, Names))

key <- c("CHROM","pos","sample")

df <- orig %>%
  select(all_of(key), Groups_o=Groups, Haps_o=Haps, Err_o=Err, Names_o=Names) %>%
  inner_join(resh %>% select(all_of(key), Groups_r=Groups, Haps_r=Haps, Err_r=Err, Names_r=Names), by=key)

check_row <- function(Names_o, Err_o, Haps_o, Names_r, Err_r, Haps_r) {
  # Assume both Err are matrices and Haps are numeric named vectors; compare after aligning by Names_o
  ok_err <- FALSE; ok_haps <- FALSE; ok_names <- FALSE

  if (is.character(Names_o) && is.character(Names_r)) {
    ok_names <- identical(Names_o, Names_r)
  }

  if (is.matrix(Err_o) && is.matrix(Err_r)) {
    Err_ra <- Err_r[Names_o, Names_o, drop=FALSE]
    ok_err <- max(abs(Err_ra - Err_o), na.rm = TRUE) <= tol
  }

  if (is.numeric(Haps_o) && is.numeric(Haps_r)) {
    Haps_ra <- Haps_r[Names_o]
    ok_haps <- max(abs(Haps_ra - Haps_o), na.rm = TRUE) <= tol
  }

  tibble(err_equal = ok_err, haps_equal = ok_haps, names_equal = ok_names)
}

res <- purrr::pmap_dfr(df, ~check_row(..5, ..3, ..4, ..8, ..6, ..7))
out <- bind_cols(df %>% select(all_of(key)), res)

cat("Total rows:", nrow(out), "\n")
cat("Err equal:", sum(out$err_equal), "\n")
cat("Haps equal:", sum(out$haps_equal), "\n")
cat("Names equal:", sum(out$names_equal), "\n")
cat("All equal:", sum(out$err_equal & out$haps_equal & out$names_equal), "\n")

mis <- out %>% filter(!(err_equal & haps_equal & names_equal))
csv_file <- paste0("mismatches_strict_", chr, ".csv")
readr::write_csv(mis, csv_file)
cat("Wrote mismatches to:", csv_file, " (", nrow(mis), " rows)\n", sep="")


