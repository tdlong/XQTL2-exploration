#!/usr/bin/env Rscript

# Compare production result vs direct est_haps_var run for one position/sample
# Usage: Rscript scripts/ErrMatrix/compare_single_position_against_prod.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(digest)
})

chr <- "chr3R"
position <- 19610000
sample_name <- "Rep01_W_F"
prod_file <- file.path("process", "ZINC2", "haplotype_results_list_format", "adapt_h4", paste0("R.haps.", chr, ".out.rds"))
test_payload <- paste0("position_data_", chr, "_", position, "_", sample_name, ".RDS")

cat("=== SINGLE-POSITION PRODUCTION VS TEST COMPARISON ===\n")
cat("Production file:", prod_file, "\n")
cat("Test payload:", test_payload, "\n\n")

stopifnot(file.exists(prod_file))
stopifnot(file.exists(test_payload))

prod_data <- readRDS(prod_file)

prod_entry <- prod_data %>%
  dplyr::filter(CHROM == chr, pos == position) %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names)) %>%
  dplyr::filter(sample == sample_name)

if (nrow(prod_entry) == 0) stop("Position/sample not found in production output")

prod_groups <- prod_entry$Groups[[1]]
prod_haps <- prod_entry$Haps[[1]]
prod_err <- prod_entry$Err[[1]]
prod_names <- prod_entry$Names[[1]]

cat("Production: Groups:", paste(prod_groups, collapse=","), " Names:", paste(prod_names, collapse=","), "\n")
cat("Production: trace(Err)=", sum(diag(prod_err), na.rm=TRUE), "\n\n")

# Run probe (ensures same code path and logs digests)
cmd <- paste(
  "Rscript",
  shQuote(file.path("scripts", "ErrMatrix", "est_haps_var_probe.R")),
  shQuote(test_payload),
  shQuote(sample_name)
)
cat("Running probe:\n", cmd, "\n\n")
status <- system(cmd)
if (status != 0) stop("Probe failed")

# Read probe report
report_file <- sub("\\.RDS$", "__probe.txt", test_payload)
cat("Reading probe report:", report_file, "\n")
stopifnot(file.exists(report_file))
report <- readLines(report_file)
cat(paste(report, collapse="\n"), "\n\n")

# Also run est_haps_var here to pull the object for matrix compare
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

payload <- readRDS(test_payload)
test_res <- est_haps_var(
  testing_position = payload$testing_position,
  sample_name = sample_name,
  df3 = payload$df4,
  founders = payload$founders,
  h_cutoff = payload$h_cutoff,
  method = payload$method,
  chr = payload$chr,
  verbose = 0
)

test_groups <- test_res$Groups
test_names <- test_res$Names
test_err <- test_res$Err

cat("Test: Groups:", paste(test_groups, collapse=","), " Names:", paste(test_names, collapse=","), "\n")
cat("Test: trace(Err)=", sum(diag(test_err), na.rm=TRUE), "\n\n")

cat("=== COMPARISON ===\n")
cat("Groups match:", identical(sort(test_groups), sort(prod_groups)), "\n")
cat("Names match:", identical(sort(test_names), sort(prod_names)), "\n")

if (is.matrix(test_err) && is.matrix(prod_err)) {
  # Align
  test_err_aligned <- test_err[prod_names, prod_names]
  diff_mat <- test_err_aligned - prod_err
  cat("Err max abs diff:", max(abs(diff_mat), na.rm=TRUE), "\n")
  cat("Err sum abs diff:", sum(abs(diff_mat), na.rm=TRUE), "\n")
} else {
  cat("One of the Err objects is not a matrix.\n")
}

cat("\nDone.\n")


