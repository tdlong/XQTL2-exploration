#!/usr/bin/env Rscript

library(tidyverse)

# Define file paths
list_results_dir <- "process/ZINC2/haplotype_results_list_format"
orig_adaptive_file <- file.path(list_results_dir, "adaptive_window_h4_results_chr3R.RDS")
smooth_original_file <- file.path(list_results_dir, "smooth_h4_results_chr3R.RDS")
reshaped_smooth_file <- file.path(list_results_dir, "smooth_h4", "R.haps.chr3R.out.rds")
reshaped_adaptive_file <- file.path(list_results_dir, "adapt_h4", "R.haps.chr3R.out.rds")

# Load files
orig_adaptive <- readRDS(orig_adaptive_file)
smooth_original <- readRDS(smooth_original_file)
reshaped_smooth <- readRDS(reshaped_smooth_file)
reshaped_adaptive <- readRDS(reshaped_adaptive_file)

# Get unique positions from each file
orig_pos <- sort(unique(orig_adaptive$pos))
smooth_pos <- sort(unique(smooth_original$pos))
reshaped_smooth_pos <- sort(unique(reshaped_smooth$pos))
reshaped_adaptive_pos <- sort(unique(reshaped_adaptive$pos))

cat("POSITION COUNTS:\n")
cat("Original adaptive:", length(orig_pos), "\n")
cat("Smooth original:", length(smooth_pos), "\n")
cat("Reshaped smooth:", length(reshaped_smooth_pos), "\n")
cat("Reshaped adaptive:", length(reshaped_adaptive_pos), "\n\n")

# Find missing positions
missing_from_smooth <- setdiff(orig_pos, smooth_pos)
missing_from_reshaped_smooth <- setdiff(orig_pos, reshaped_smooth_pos)
missing_from_reshaped_adaptive <- setdiff(orig_pos, reshaped_adaptive_pos)

cat("MISSING FROM SMOOTH ORIGINAL (", length(missing_from_smooth), " positions):\n")
if (length(missing_from_smooth) > 0) {
  cat(paste(missing_from_smooth, collapse=", "), "\n")
}
cat("\n")

cat("MISSING FROM RESHAPED SMOOTH (", length(missing_from_reshaped_smooth), " positions):\n")
if (length(missing_from_reshaped_smooth) > 0) {
  cat(paste(missing_from_reshaped_smooth, collapse=", "), "\n")
}
cat("\n")

cat("MISSING FROM RESHAPED ADAPTIVE (", length(missing_from_reshaped_adaptive), " positions):\n")
if (length(missing_from_reshaped_adaptive) > 0) {
  cat(paste(missing_from_reshaped_adaptive, collapse=", "), "\n")
}
cat("\n")

# Show first and last few positions of each
cat("FIRST 10 POSITIONS IN EACH FILE:\n")
cat("Original adaptive:", paste(head(orig_pos, 10), collapse=", "), "\n")
cat("Smooth original:", paste(head(smooth_pos, 10), collapse=", "), "\n")
cat("Reshaped smooth:", paste(head(reshaped_smooth_pos, 10), collapse=", "), "\n")
cat("Reshaped adaptive:", paste(head(reshaped_adaptive_pos, 10), collapse=", "), "\n\n")

cat("LAST 10 POSITIONS IN EACH FILE:\n")
cat("Original adaptive:", paste(tail(orig_pos, 10), collapse=", "), "\n")
cat("Smooth original:", paste(tail(smooth_pos, 10), collapse=", "), "\n")
cat("Reshaped smooth:", paste(tail(reshaped_smooth_pos, 10), collapse=", "), "\n")
cat("Reshaped adaptive:", paste(tail(reshaped_adaptive_pos, 10), collapse=", "), "\n")
