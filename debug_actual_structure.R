#!/usr/bin/env Rscript

library(tidyverse)

# Load files
orig <- readRDS('haplotype_results_list_format/adaptive_window_h4_results_chr3R.RDS')
resh <- readRDS('haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds')

cat("=== ORIGINAL FILE ===\n")
cat("Class:", class(orig), "\n")
cat("Dimensions:", nrow(orig), "x", ncol(orig), "\n")
cat("Columns:", paste(names(orig), collapse=", "), "\n")
cat("First row sample class:", class(orig$sample[1]), "length:", length(orig$sample[1]), "\n")
cat("First row Names class:", class(orig$Names[1]), "length:", length(orig$Names[1]), "\n")
cat("First row Err class:", class(orig$Err[1]), "length:", length(orig$Err[1]), "\n")
cat("First Err[[1]] class:", class(orig$Err[[1]]), "is.matrix:", is.matrix(orig$Err[[1]]), "\n")
if (is.matrix(orig$Err[[1]])) {
  cat("First Err[[1]] dim:", paste(dim(orig$Err[[1]]), collapse="x"), "\n")
  cat("First Err[[1]] rownames:", paste(rownames(orig$Err[[1]]), collapse=","), "\n")
}

cat("\n=== RESHAPED FILE ===\n")
cat("Class:", class(resh), "\n")
cat("Dimensions:", nrow(resh), "x", ncol(resh), "\n")
cat("Columns:", paste(names(resh), collapse=", "), "\n")
cat("First row sample class:", class(resh$sample[1]), "length:", length(resh$sample[1]), "\n")
cat("First row Names class:", class(resh$Names[1]), "length:", length(resh$Names[1]), "\n")
cat("First row Err class:", class(resh$Err[1]), "length:", length(resh$Err[1]), "\n")
cat("First Err[[1]] class:", class(resh$Err[[1]]), "is.matrix:", is.matrix(resh$Err[[1]]), "\n")
if (is.matrix(resh$Err[[1]])) {
  cat("First Err[[1]] dim:", paste(dim(resh$Err[[1]]), collapse="x"), "\n")
  cat("First Err[[1]] rownames:", paste(rownames(resh$Err[[1]]), collapse=","), "\n")
}
