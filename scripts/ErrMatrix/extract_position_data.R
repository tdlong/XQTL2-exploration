#!/usr/bin/env Rscript

# Extract data for a specific position (df4) that gets passed to the haplotype estimator
# Usage: Rscript scripts/ErrMatrix/extract_position_data.R <chr> <position> <output_dir> <param_file>

library(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript scripts/ErrMatrix/extract_position_data.R <chr> <position> <output_dir> <param_file>\n")
  cat("Example: Rscript scripts/ErrMatrix/extract_position_data.R chr3R 19610000 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R\n")
  quit(status = 1)
}

chr <- args[1]
position <- as.numeric(args[2])
output_dir <- args[3]
param_file <- args[4]

cat("=== EXTRACTING POSITION DATA ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", position, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Load parameters
source(param_file, local = TRUE)
founders <- get("founders")

cat("Founders:", paste(founders, collapse=", "), "\n\n")

# Load RefAlt data
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
cat("Loading RefAlt data from:", refalt_file, "\n")

if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

# Load data using the same approach as BASE_VAR_WIDE.R
df <- read.table(refalt_file, header = TRUE)
cat("Loaded", nrow(df), "rows from RefAlt file\n")

# Transform to long format (same as BASE_VAR_WIDE.R)
df2 <- df %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    N = REF + ALT
  ) %>%
  select(-REF, -ALT) %>%
  as_tibble()

# Apply quality filter (founders must be fixed)
founder_wide <- df2 %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

quality_filtered_positions <- founder_wide %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  ) %>%
  pull(POS)

# Filter to quality positions
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

cat("Quality-filtered positions:", length(quality_filtered_positions), "\n")
cat("Processed", nrow(df3), "rows for", length(unique(df3$name)), "samples\n\n")

# Create subset for the target position (500kb window)
window_size <- 500000
df4 <- df3 %>%
  filter(POS >= position - window_size/2 & 
         POS <= position + window_size/2)

cat("Subsetted data for position", position, ":\n")
cat("  Window size:", window_size, "bp\n")
cat("  SNPs in window:", nrow(df4), "\n")
cat("  Columns:", paste(names(df4), collapse=", "), "\n")
cat("  Position range:", min(df4$POS), "to", max(df4$POS), "\n\n")

# Save the subsetted data
output_file <- paste0("position_data_", chr, "_", position, ".RDS")
saveRDS(df4, output_file)
cat("Saved subsetted data to:", output_file, "\n")

# Also save some metadata
metadata <- list(
  chr = chr,
  position = position,
  window_size = window_size,
  n_snps = nrow(df4),
  founders = founders,
  samples = unique(df4$name),
  pos_range = c(min(df4$POS), max(df4$POS))
)

metadata_file <- paste0("position_metadata_", chr, "_", position, ".RDS")
saveRDS(metadata, metadata_file)
cat("Saved metadata to:", metadata_file, "\n")

cat("\nDone.\n")
