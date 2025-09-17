#!/usr/bin/env Rscript

# Extract Testing Positions for Bug Analysis
# This script extracts the 7 testing positions around 19780000-19840000 on chr3R
# from both adaptive and fixed window methods for comparison

# Load required libraries
library(tidyverse)

# Define the testing positions (every 10000 on the 10000)
testing_positions <- c(19780000, 19790000, 19800000, 19810000, 19820000, 19830000, 19840000)

# File paths on the server
adapt_file <- "/dfs7/adl/tdlong/exploration/XQTL2-exploration/process/ZINC2/haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds"
fixed_file <- "/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/R.haps.chr3R.out.rds"

# Output file path (in current working directory)
output_file <- "testing_positions_comparison.rds"

cat("Extracting testing positions for bug analysis...\n")
cat("Testing positions:", paste(testing_positions, collapse = ", "), "\n\n")

# Function to extract positions from a file
extract_positions <- function(file_path, method_name) {
  cat("Loading", method_name, "data from:", basename(file_path), "\n")
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Load the data
  data <- readRDS(file_path)
  
  # Check if data is a list or data frame
  if (is.list(data) && !is.data.frame(data)) {
    # If it's a list, we need to find the data frame
    # Look for common column names that might indicate the main data
    possible_data <- NULL
    for (i in seq_along(data)) {
      if (is.data.frame(data[[i]]) && "pos" %in% names(data[[i]])) {
        possible_data <- data[[i]]
        break
      }
    }
    if (is.null(possible_data)) {
      stop("Could not find data frame with 'pos' column in list structure")
    }
    data <- possible_data
  }
  
  # Check if 'pos' column exists
  if (!"pos" %in% names(data)) {
    stop("No 'pos' column found in data. Available columns: ", paste(names(data), collapse = ", "))
  }
  
  # Extract the testing positions
  extracted <- data %>%
    dplyr::filter(pos %in% testing_positions) %>%
    dplyr::arrange(pos) %>%
    dplyr::mutate(method = method_name)
  
  cat("Found", nrow(extracted), "positions out of", length(testing_positions), "requested\n")
  cat("Positions found:", paste(extracted$pos, collapse = ", "), "\n")
  
  return(extracted)
}

# Extract data from both methods
cat(strrep("=", 60), "\n")
cat("EXTRACTING ADAPTIVE METHOD DATA\n")
cat(strrep("=", 60), "\n")
adapt_data <- extract_positions(adapt_file, "adapt")

cat("\n")
cat(strrep("=", 60), "\n")
cat("EXTRACTING FIXED METHOD DATA\n")
cat(strrep("=", 60), "\n")
fixed_data <- extract_positions(fixed_file, "fixed")

# Combine the data
cat("\n")
cat(strrep("=", 60), "\n")
cat("COMBINING DATA\n")
cat(strrep("=", 60), "\n")

# Check if both datasets have the same columns (excluding method)
adapt_cols <- names(adapt_data)[names(adapt_data) != "method"]
fixed_cols <- names(fixed_data)[names(fixed_data) != "method"]

if (!setequal(adapt_cols, fixed_cols)) {
  cat("WARNING: Column names differ between methods\n")
  cat("Adaptive columns:", paste(adapt_cols, collapse = ", "), "\n")
  cat("Fixed columns:", paste(fixed_cols, collapse = ", "), "\n")
  cat("Common columns:", paste(intersect(adapt_cols, fixed_cols), collapse = ", "), "\n")
  
  # Use only common columns
  common_cols <- intersect(adapt_cols, fixed_cols)
  adapt_data <- adapt_data[, c(common_cols, "method")]
  fixed_data <- fixed_data[, c(common_cols, "method")]
}

# Combine the data
combined_data <- dplyr::bind_rows(adapt_data, fixed_data) %>%
  dplyr::arrange(pos, method)

cat("Combined data shape:", nrow(combined_data), "rows,", ncol(combined_data), "columns\n")
cat("Methods:", paste(unique(combined_data$method), collapse = ", "), "\n")
cat("Positions:", paste(unique(combined_data$pos), collapse = ", "), "\n")

# Save the combined data
cat("\n")
cat(strrep("=", 60), "\n")
cat("SAVING RESULTS\n")
cat(strrep("=", 60), "\n")

saveRDS(combined_data, output_file)
cat("Saved combined data to:", output_file, "\n")

# Display summary
cat("\n")
cat(strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")

# Show data structure
cat("Data structure:\n")
str(combined_data)

# Show first few rows
cat("\nFirst few rows:\n")
print(head(combined_data, 10))

# Show position counts by method
cat("\nPosition counts by method:\n")
position_counts <- combined_data %>%
  dplyr::count(method, pos) %>%
  dplyr::arrange(method, pos)
print(position_counts)

cat("\nExtraction complete! File saved as:", output_file, "\n")
cat("You can now download this file to your local machine for analysis.\n")
