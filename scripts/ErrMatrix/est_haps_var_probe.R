#!/usr/bin/env Rscript

# Probe wrapper to call production est_haps_var and emit strong verification diagnostics.
# Usage:
#   Rscript scripts/ErrMatrix/est_haps_var_probe.R <data_payload_rds> <sample_name>
# Where <data_payload_rds> is the file produced by extract_single_position_data.R

suppressPackageStartupMessages({
  library(digest)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript scripts/ErrMatrix/est_haps_var_probe.R <data_payload_rds> <sample_name>\n")
  quit(status = 1)
}

data_file <- args[1]
sample_name <- args[2]

payload <- readRDS(data_file)
df4 <- payload$df4
testing_position <- payload$testing_position
founders <- payload$founders
h_cutoff <- payload$h_cutoff
method <- payload$method
chr <- payload$chr

# Hash inputs
hash_df4 <- digest(df4, algo = "sha1")
hash_founders <- digest(founders, algo = "sha1")
hash_params <- digest(list(testing_position, sample_name, founders, h_cutoff, method, chr), algo = "sha1")

cat("=== PROBE: INPUT DIGESTS ===\n")
cat("df4_sha1:", hash_df4, "\n")
cat("founders_sha1:", hash_founders, "\n")
cat("params_sha1:", hash_params, "\n\n")

# Source production code (do not change production file)
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Hash code for function identity
code_hash <- tryCatch(digest(deparse(body(est_haps_var)), algo = "sha1"), error = function(e) NA)
cat("est_haps_var_code_sha1:", code_hash, "\n\n")

cat("Running est_haps_var...\n")
result <- est_haps_var(
  testing_position = testing_position,
  sample_name = sample_name,
  df3 = df4,
  founders = founders,
  h_cutoff = h_cutoff,
  method = method,
  chr = chr,
  verbose = 1
)

Groups <- result$Groups
Haps <- result$Haps
Err <- result$Err
Names <- result$Names

trace_Err <- tryCatch(sum(diag(Err), na.rm = TRUE), error = function(e) NA)
fro_Err <- tryCatch(sqrt(sum(Err^2, na.rm = TRUE)), error = function(e) NA)

cat("\n=== RESULT DIAGNOSTICS ===\n")
cat("Groups:", paste(Groups, collapse=","), "\n")
cat("Names:", paste(Names, collapse=","), "\n")
cat("trace(Err):", format(trace_Err, digits = 6), "\n")
cat("||Err||_F:", format(fro_Err, digits = 6), "\n")

# Save a compact report alongside data file
out_txt <- sub("\\.RDS$", "__probe.txt", data_file)
sink(out_txt)
cat("df4_sha1:", hash_df4, "\n")
cat("founders_sha1:", hash_founders, "\n")
cat("params_sha1:", hash_params, "\n")
cat("est_haps_var_code_sha1:", code_hash, "\n")
cat("Groups:", paste(Groups, collapse=","), "\n")
cat("Names:", paste(Names, collapse=","), "\n")
cat("trace_Err:", trace_Err, "\n")
cat("fro_Err:", fro_Err, "\n")
sink()
cat("\nProbe report saved:", out_txt, "\n")


