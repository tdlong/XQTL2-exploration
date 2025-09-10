#!/usr/bin/env Rscript

# Convert haplotype estimation results to list-based format
# This script transforms the current wide format to the desired list format
# where each position has one row with list columns for samples, groups, haps, err, and names

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript convert_haplotype_to_list_format.R <chr> <method> <parameter> <output_dir>")
}

chr <- args[1]
method <- args[2]
parameter <- args[3]
output_dir <- args[4]

cat("=== CONVERTING HAPLOTYPE TO LIST FORMAT ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n\n")

# Determine the input file based on method and parameter
results_dir <- file.path(output_dir, "haplotype_results")

if (method == "fixed") {
  input_file <- file.path(results_dir, paste0("fixed_window_", parameter, "kb_results_", chr, ".RDS"))
  output_file <- file.path(results_dir, paste0("fixed_window_", parameter, "kb_list_format_", chr, ".RDS"))
} else if (method == "adaptive") {
  input_file <- file.path(results_dir, paste0("adaptive_window_h", parameter, "_results_", chr, ".RDS"))
  output_file <- file.path(results_dir, paste0("adaptive_window_h", parameter, "_list_format_", chr, ".RDS"))
} else if (method == "smooth_h4") {
  input_file <- file.path(results_dir, paste0("smooth_h4_results_", chr, ".RDS"))
  output_file <- file.path(results_dir, paste0("smooth_h4_list_format_", chr, ".RDS"))
} else {
  stop("Unknown method: ", method)
}

# Check if input file exists
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

cat("✓ Loading input file:", basename(input_file), "\n")
haplotype_data <- readRDS(input_file)

cat("Input data dimensions:", nrow(haplotype_data), "rows,", ncol(haplotype_data), "columns\n")
cat("Columns:", paste(names(haplotype_data), collapse = ", "), "\n\n")

# Get founder names (B1, B2, B3, B4, B5, B6, B7, AB8)
founder_cols <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
founder_names <- founder_cols[founder_cols %in% names(haplotype_data)]

if (length(founder_names) == 0) {
  stop("No founder columns found. Expected columns: B1, B2, B3, B4, B5, B6, B7, AB8")
}

cat("Found founder columns:", paste(founder_names, collapse = ", "), "\n")

# Get unique samples
samples <- unique(haplotype_data$sample)
cat("Samples:", length(samples), "\n")
cat("Sample names:", paste(samples, collapse = ", "), "\n\n")

# Get unique positions
positions <- unique(haplotype_data$pos)
cat("Unique positions:", length(positions), "\n")
cat("Position range:", min(positions), "-", max(positions), "\n\n")

# Function to create groups (for now, each founder gets its own group)
create_groups <- function(founder_names) {
  groups <- 1:length(founder_names)
  names(groups) <- founder_names
  return(groups)
}

# Function to create error matrix (placeholder - we don't have the actual error matrix)
create_error_matrix <- function(founder_names) {
  n_founders <- length(founder_names)
  # Create a simple diagonal matrix as placeholder
  error_matrix <- diag(0.01, n_founders)
  rownames(error_matrix) <- founder_names
  colnames(error_matrix) <- founder_names
  return(error_matrix)
}

# Convert to list format
cat("Converting to list format...\n")

list_format_data <- positions %>%
  map_dfr(function(pos) {
    # Get data for this position
    pos_data <- haplotype_data %>%
      filter(pos == !!pos) %>%
      arrange(sample)
    
    if (nrow(pos_data) == 0) {
      return(NULL)
    }
    
    # Create groups (same for all samples at this position)
    groups <- create_groups(founder_names)
    
    # Create error matrix (placeholder)
    error_matrix <- create_error_matrix(founder_names)
    
    # Create list of samples
    sample_list <- list(pos_data$sample)
    
    # Create list of groups (same for all samples)
    groups_list <- list(rep(list(groups), length(samples)))
    
    # Create list of haplotype frequencies for each sample
    haps_list <- list()
    for (i in 1:nrow(pos_data)) {
      sample_haps <- pos_data[i, founder_names] %>%
        as.numeric() %>%
        set_names(founder_names)
      haps_list[[i]] <- sample_haps
    }
    haps_list <- list(haps_list)
    
    # Create list of error matrices (same for all samples)
    err_list <- list(rep(list(error_matrix), length(samples)))
    
    # Create list of names (same for all samples)
    names_list <- list(rep(list(founder_names), length(samples)))
    
    # Return as tibble row
    tibble(
      CHROM = chr,
      pos = pos,
      sample = sample_list,
      Groups = groups_list,
      Haps = haps_list,
      Err = err_list,
      Names = names_list
    )
  })

cat("✓ Conversion complete!\n")
cat("Output data dimensions:", nrow(list_format_data), "rows,", ncol(list_format_data), "columns\n")
cat("Columns:", paste(names(list_format_data), collapse = ", "), "\n\n")

# Show sample of the output
cat("Sample of converted data:\n")
print(head(list_format_data, 3))

# Show structure of list columns
if (nrow(list_format_data) > 0) {
  cat("\nStructure of list columns:\n")
  cat("Sample list length:", length(list_format_data$sample[[1]]), "\n")
  cat("Groups list length:", length(list_format_data$Groups[[1]]), "\n")
  cat("Haps list length:", length(list_format_data$Haps[[1]]), "\n")
  cat("Err list length:", length(list_format_data$Err[[1]]), "\n")
  cat("Names list length:", length(list_format_data$Names[[1]]), "\n")
  
  # Show example of first position
  if (nrow(list_format_data) > 0) {
    cat("\nExample - First position:\n")
    cat("Samples:", paste(list_format_data$sample[[1]], collapse = ", "), "\n")
    cat("Groups:", paste(list_format_data$Groups[[1]][[1]], collapse = ", "), "\n")
    cat("Haps (first sample):", paste(round(list_format_data$Haps[[1]][[1]], 4), collapse = ", "), "\n")
    cat("Names:", paste(list_format_data$Names[[1]][[1]], collapse = ", "), "\n")
  }
}

# Save the converted data
cat("\nSaving converted data to:", basename(output_file), "\n")
saveRDS(list_format_data, output_file)

cat("✓ File saved successfully!\n")
cat("File size:", round(file.size(output_file) / 1024^2, 2), "MB\n")

cat("\n=== CONVERSION COMPLETE ===\n")
cat("The converted file is ready for use with your other analysis pipeline.\n")
cat("Note: Error matrices are currently placeholders (diagonal matrices).\n")
cat("If you need the actual error matrices, we'll need to modify the haplotype estimation process.\n")
