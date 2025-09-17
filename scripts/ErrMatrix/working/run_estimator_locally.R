#!/usr/bin/env Rscript

# Local wrapper: load df3 + args and call est_haps_var() exactly

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/ErrMatrix/run_estimator_locally.R <payload_rds> <base_var_wide_path>")
}

payload_rds <- args[1]
base_var_wide_path <- args[2]

if (!file.exists(payload_rds)) stop("Payload not found: ", payload_rds)
if (!file.exists(base_var_wide_path)) stop("BASE_VAR_WIDE.R not found: ", base_var_wide_path)

payload <- readRDS(payload_rds)

# Bring in est_haps_var() definition
source(base_var_wide_path)

df3 <- payload$df3
est_args <- payload$args

res <- do.call(est_haps_var, c(list(df3 = df3), est_args))

print(str(res))

saveRDS(res, sub("\\.rds$", ".result.rds", payload_rds))
cat("Saved:", sub("\\.rds$", ".result.rds", payload_rds), "\n")


