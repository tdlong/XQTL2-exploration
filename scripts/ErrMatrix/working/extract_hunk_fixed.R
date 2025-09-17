#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# Include process_refalt_data function so script is self-contained - EXACT COPY FROM PRODUCTION
process_refalt_data <- function(refalt_file, founders) {
  # Load RefAlt data and process into df3 format - EXACT from working code
  cat("Loading RefAlt data from:", refalt_file, "\n")
  
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Process into df3 format (one row per sample per position) - EXACT from working code
  df3 <- refalt_data %>%
    pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
    mutate(
      RefAlt = str_sub(lab, 1, 3),
      name = str_sub(lab, 5)
    ) %>%
    dplyr::select(-lab) %>%
    { cat("After pivot_longer, columns:", paste(names(.), collapse=", "), "\n"); . } %>%
    pivot_wider(names_from = RefAlt, values_from = count) %>%
    { cat("After pivot_wider, columns:", paste(names(.), collapse=", "), "\n"); . } %>%
    mutate(
      freq = REF / (REF + ALT),
      N = REF + ALT
    ) %>%
    # Debug: check what columns exist
    { cat("Columns before select:", paste(names(.), collapse=", "), "\n"); . } %>%
    dplyr::select(-REF, -ALT)
  
  # Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
  founder_wide <- df3 %>%
    filter(name %in% founders) %>%
    dplyr::select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  quality_filtered_positions <- founder_wide %>%
    filter(
      if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
    ) %>%
    pull(POS)
  
  # Filter to quality positions and include sample data
  df3 <- df3 %>%
    filter(POS %in% quality_filtered_positions)
  
  # Convert to wide format before returning - ONE LINER OPTIMIZATION
  df3 <- df3 %>% dplyr::select(POS, name, freq) %>% pivot_wider(names_from = name, values_from = freq)
  
  cat("✓ Processed", nrow(df3), "rows for", ncol(df3)-1, "samples in WIDE format\n")
  return(df3)
}

extract_hunk <- function(chr, test_position, test_sample, method="adaptive", parameter=4, param_file, output_dir="process/ZINC2") {
  cat("=== EXTRACTING DATA FOR SINGLE POSITION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Position:", test_position, "\n")
  cat("Sample:", test_sample, "\n")
  cat("Method:", method, "\n")
  cat("Parameter:", parameter, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n\n")
  
  # Load parameters (defines founders and other variables)
  source(param_file, local = TRUE)
  
  # Find the RefAlt file for this chromosome
  refalt_file <- file.path(output_dir, paste0("REFALT.", chr, ".txt"))
  
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  cat("Using RefAlt file:", refalt_file, "\n")
  cat("Founders:", paste(founders, collapse = ", "), "\n\n")
  
  # Process RefAlt data into df3 format (EXACT same as production)
  df3 <- process_refalt_data(refalt_file, founders)
  
  # Filter to the specific position and sample
  df4 <- df3 %>%
    filter(POS == test_position) %>%
    dplyr::select(-POS) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "freq") %>%
    filter(sample == test_sample)
  
  if (nrow(df4) == 0) {
    stop("Sample ", test_sample, " not found at position ", test_position)
  }
  
  cat("Found sample", test_sample, "at position", test_position, "\n")
  cat("df3 dimensions:", nrow(df3), "x", ncol(df3), "\n")
  cat("df4 dimensions:", nrow(df4), "x", ncol(df4), "\n\n")
  
  # Package the data and arguments for est_haps_var
  payload <- list(
    df3 = df3,  # The full df3 for the position
    args = list(
      testing_position = test_position,
      sample_name = test_sample,
      founders = founders,
      h_cutoff = parameter,
      method = method,
      window_size_bp = 100000,  # Default window size
      chr = chr,
      verbose = 3  # Maximum verbosity for debugging
    )
  )
  
  # Save the payload
  output_file <- paste0("hunk_data_", chr, "_", test_position, "_", test_sample, ".rds")
  saveRDS(payload, output_file)
  
  cat("✓ Saved hunk data to:", output_file, "\n")
  cat("✓ Ready for debugging with exact production data\n")
}

# Run if called directly
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 4) {
    chr <- args[1]
    test_position <- as.numeric(args[2])
    test_sample <- args[3]
    parameter <- as.numeric(args[4])
    param_file <- if (length(args) > 4) args[5] else "helpfiles/ZINC2_haplotype_parameters.R"
    output_dir <- if (length(args) > 5) args[6] else "process/ZINC2"
    
    extract_hunk(chr, test_position, test_sample, parameter=parameter, param_file=param_file, output_dir=output_dir)
  } else {
    cat("Usage: Rscript extract_hunk.r <chr> <position> <sample> <h_cutoff> [param_file] [output_dir]\n")
    cat("Example: Rscript extract_hunk.r chr3R 19780000 Rep01_W_F 4\n")
  }
}
