#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

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
  source(param_file)
  if (!exists("founders")) {
    stop("Param file must define 'founders' variable")
  }
  cat("Founders:", paste(founders, collapse = ", "), "\n")
  
  # Load RefAlt data using production function
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  df3 <- process_refalt_data(refalt_file, founders)
  cat("✓ df3 loaded:", nrow(df3), "rows x", ncol(df3), "cols\n")
  
  # Pre-subset df3 to df4 for this position (testing_position ± 500kb)
  max_window <- 500000  # Largest window size in est_haps_var
  window_start <- max(1, test_position - max_window/2)
  window_end <- test_position + max_window/2
      
  df4 <- df3 %>% 
    dplyr::filter(POS >= window_start & POS <= window_end)
  
  cat("✓ df4 subsetted:", nrow(df4), "SNPs in window", window_start, "-", window_end, "\n")
  
  # Check if sample exists in the data
  available_samples <- names(df4)[!names(df4) %in% c("POS")]
  if (!test_sample %in% available_samples) {
    stop("Sample '", test_sample, "' not found in data. Available samples: ", 
         paste(head(available_samples, 10), collapse = ", "), 
         if(length(available_samples) > 10) "..." else "")
  }
  cat("✓ Sample", test_sample, "found in data\n")
  
  # Package data for debugging
  payload <- list(
    df3 = df4,  # Pre-subsetted data like production
    args = list(
      testing_position = test_position,
      sample_name = test_sample,
      founders = founders,
      h_cutoff = parameter,
      method = method,
      window_size_bp = NULL,
      chr = chr,
      verbose = 2
    )
  )
  
  outfile <- sprintf("hunk_data_%s_%d_%s_h%d.rds", chr, test_position, test_sample, parameter)
  saveRDS(payload, outfile)
  cat("Saved:", outfile, "\n")
  cat("=== Extraction complete ===\n")
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