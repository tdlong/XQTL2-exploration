#!/usr/bin/env Rscript

# Process centromere SNP files for XQTL2 exploration
# This script reads the 2LRhet and 3LRhet files and creates a combined dataframe

library(tidyverse)

# Function to process a single file
process_centromere_file <- function(file_path) {
  # Read the second line (skip first line with count)
  snp_line <- readLines(file_path, n = 2)[2]
  
  # Split by whitespace and remove empty strings
  snp_positions <- strsplit(snp_line, "\\s+")[[1]]
  snp_positions <- snp_positions[snp_positions != ""]
  
  # Remove the first element (header like "AA_header1")
  snp_positions <- snp_positions[-1]
  
  # Split each position on colon
  snp_data <- strsplit(snp_positions, ":")
  
  # Create dataframe
  df <- data.frame(
    CHROM = sapply(snp_data, function(x) paste0("chr", x[1])),
    pos = as.numeric(sapply(snp_data, function(x) x[2])),
    stringsAsFactors = FALSE
  )
  
  return(df)
}

# Process both files
cat("Processing 2LRhet file...\n")
df_2lr <- process_centromere_file("../helpfiles/2LRhet_ZIinbred_ADLs_I80_Q50_mm1_mac3_snps_mskn40_hap.txt")

cat("Processing 3LRhet file...\n")
df_3lr <- process_centromere_file("../helpfiles/3LRhet_ZIinbred_ADLs_I80_Q50_mm1_mac3_snps_mskn40_hap.txt")

# Combine the dataframes
cat("Combining dataframes...\n")
combined_df <- rbind(df_2lr, df_3lr)

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Total positions:", nrow(combined_df), "\n")
cat("Chromosome levels:\n")
print(table(combined_df$CHROM))

# Detailed summary with min/max positions and range
summary_df <- combined_df %>%
  group_by(CHROM) %>%
  summarise(
    n_positions = n(),
    min_pos = min(pos),
    max_pos = max(pos),
    range_mb = (max(pos) - min(pos)) / 1e6,
    .groups = "drop"
  )

cat("\nPositions per chromosome:\n")
print(summary_df)

# Save the result
output_file <- "sasha_good.rds"
cat("\nSaving to", output_file, "...\n")
saveRDS(combined_df, output_file)

cat("Done! Saved", nrow(combined_df), "SNP positions to", output_file, "\n")
