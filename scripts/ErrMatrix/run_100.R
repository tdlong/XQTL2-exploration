#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

source("scripts/ErrMatrix/haplotype_error_workbench.R")

args <- commandArgs(trailingOnly = TRUE)
h_cutoff <- if (length(args) >= 1) as.numeric(args[[1]]) else 4

run_batch_100_with_summary(h_cutoff = h_cutoff)


