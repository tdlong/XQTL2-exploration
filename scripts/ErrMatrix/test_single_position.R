#!/usr/bin/env Rscript

# Test script to reproduce Friday night results for one position
# Usage: Rscript scripts/ErrMatrix/test_single_position.R

# Load the Friday night version of BASE_VAR_WIDE.R
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")

# Test parameters
chr <- "chr3R"
testing_position <- 19610000
sample_name <- "Rep01_W_F"
param_file <- "helpfiles/ZINC2_haplotype_parameters.R"
output_dir <- "process/ZINC2"

cat("=== TESTING SINGLE POSITION ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", testing_position, "\n")
cat("Sample:", sample_name, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file)
founders <- get("founders", envir = .GlobalEnv)
parameter <- get("parameter", envir = .GlobalEnv)
method <- get("method", envir = .GlobalEnv)

cat("Loaded parameters:\n")
cat("  Founders:", paste(founders, collapse=", "), "\n")
cat("  Parameter:", parameter, "\n")
cat("  Method:", method, "\n\n")

# Load RefAlt data
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
cat("Loading RefAlt data from:", refalt_file, "\n")

if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

# Process RefAlt data (this will be in wide format)
df3 <- process_refalt_data(refalt_file, founders)
cat("Processed", nrow(df3), "rows for", ncol(df3)-1, "samples\n\n")

# Create subset for this position (500kb window)
window_size <- 500000
df4 <- df3 %>%
  dplyr::filter(POS >= testing_position - window_size/2 & 
                POS <= testing_position + window_size/2)

cat("Subsetted data for position", testing_position, ":\n")
cat("  Window size:", window_size, "bp\n")
cat("  SNPs in window:", nrow(df4), "\n")
cat("  Columns:", paste(names(df4), collapse=", "), "\n\n")

# Save the subsetted data
subset_file <- paste0("test_subset_", chr, "_", testing_position, "_", sample_name, ".RDS")
saveRDS(df4, subset_file)
cat("Saved subsetted data to:", subset_file, "\n\n")

# Run the estimator
cat("Running estimator...\n")
result <- estimate_haplotypes_list_format(
  testing_position = testing_position,
  sample_name = sample_name,
  df3 = df4,
  founders = founders,
  h_cutoff = parameter,
  method = method,
  chr = chr,
  verbose = 2
)

cat("\n=== RESULT ===\n")
cat("Groups:", paste(result$Groups, collapse=","), "\n")
cat("Haplotypes:\n")
print(result$Haps)
cat("Error matrix:\n")
print(result$Err)
cat("Names:", paste(result$Names, collapse=","), "\n")

# Save the result
result_file <- paste0("test_result_", chr, "_", testing_position, "_", sample_name, ".RDS")
saveRDS(result, result_file)
cat("\nSaved result to:", result_file, "\n")

cat("\nDone.\n")
