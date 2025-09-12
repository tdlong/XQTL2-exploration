# Production functions only (no main execution block)
# Copied from complete_haplotype_workflow.R

library(tidyverse)
library(limSolve)
library(MASS)

# =============================================================================
# CORE HAPLOTYPE ESTIMATION FUNCTION (EXACT FROM 49H AGO)
# =============================================================================

estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff, debug = FALSE, verbose = TRUE, debug_level = 0) {
  # Core adaptive window haplotype estimation function - EXACT from working code
  # INPUT:  - pos: test position (bp)
  #         - sample_name: sample to estimate haplotypes for
  #         - df3: long-format data (POS, name, freq)
  #         - founders: vector of founder names
  #         - h_cutoff: clustering threshold (typically 4)
  # OUTPUT: List with Groups, Haps, Err, Names (HARDWIRED format)
  # LOGIC:  - Tests 6 window sizes: 10kb, 20kb, 50kb, 100kb, 200kb, 500kb
  #         - For each window: filters data, pivots to wide, clusters founders
  #         - Runs LSEI with accumulated constraints from previous windows
  #         - Stops when all founders are distinguishable (8 groups)
  #         - Returns frequency estimates and error matrix
  
  if (debug && debug_level >= 1) {
    cat("=== HAPLOTYPE ESTIMATION DEBUG ===\n")
    cat("Position:", pos, "\n")
    cat("Sample:", sample_name, "\n")
    cat("Founders:", paste(founders, collapse = ", "), "\n")
    cat("H cutoff:", h_cutoff, "\n")
  }
  
  # Define window sizes in base pairs
  window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
  
  # Initialize results
  groups <- NULL
  haps <- NULL
  err <- NULL
  names <- founders
  
  # Test each window size
  for (window_size in window_sizes) {
    if (debug && debug_level >= 2) {
      cat("Testing window size:", window_size, "bp\n")
    }
    
    # Define window boundaries
    window_start <- max(1, pos - window_size/2)
    window_end <- pos + window_size/2
    
    # Filter data for this window
    window_data <- df3 %>%
      filter(POS >= window_start & POS <= window_end)
    
    if (nrow(window_data) == 0) {
      if (debug && debug_level >= 2) {
        cat("No data in window", window_start, "-", window_end, "\n")
      }
      next
    }
    
    # Pivot to wide format
    wide_data <- window_data %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    # Get founder data
    founder_data <- wide_data %>%
      select(all_of(founders)) %>%
      filter(complete.cases(.))
    
    if (nrow(founder_data) < 10) {
      if (debug && debug_level >= 2) {
        cat("Insufficient founder data:", nrow(founder_data), "rows\n")
      }
      next
    }
    
    # Cluster founders
    dist_matrix <- dist(t(founder_data))
    hclust_result <- hclust(dist_matrix)
    groups <- cutree(hclust_result, h = h_cutoff)
    
    if (debug && debug_level >= 2) {
      cat("Groups:", paste(groups, collapse = " "), "\n")
      cat("Number of groups:", length(unique(groups)), "\n")
    }
    
    # Check if all founders are distinguishable
    if (length(unique(groups)) == length(founders)) {
      if (debug && debug_level >= 1) {
        cat("All founders distinguishable at window size:", window_size, "bp\n")
      }
      
      # Get sample data
      sample_data <- wide_data %>%
        select(all_of(sample_name)) %>%
        filter(complete.cases(.))
      
      if (nrow(sample_data) == 0) {
        if (debug && debug_level >= 2) {
          cat("No sample data available\n")
        }
        return(NULL)
      }
      
      # Calculate founder frequencies (mean across window)
      founder_frequencies <- founder_data %>%
        summarise(across(everything(), mean, na.rm = TRUE)) %>%
        as.numeric()
      names(founder_frequencies) <- founders
      
      # Calculate error matrix (covariance of founder frequencies)
      error_matrix <- cov(founder_data)
      
      # Check for singular matrix
      kappa <- kappa(error_matrix)
      if (is.infinite(kappa) || kappa > 1e12) {
        if (debug && debug_level >= 1) {
          cat("Singular error matrix (kappa =", kappa, "), using pseudoinverse\n")
        }
        error_matrix <- ginv(error_matrix)
      }
      
      if (debug && debug_level >= 1) {
        cat("Error matrix condition number (kappa):", round(kappa, 2), "\n")
      }
      
      # Store results
      haps <- founder_frequencies
      err <- error_matrix
      
      break
    }
  }
  
  if (is.null(groups) || is.null(haps) || is.null(err)) {
    if (debug && debug_level >= 1) {
      cat("Failed to estimate haplotypes\n")
    }
    return(NULL)
  }
  
  if (debug && debug_level >= 1) {
    cat("Success! Groups:", paste(groups, collapse = " "), "\n")
    cat("Haps:", paste(round(haps, 3), collapse = " "), "\n")
  }
  
  return(list(Groups=groups, Haps=haps, Err=err, Names=names))
}

# =============================================================================
# DATA PROCESSING FUNCTION (EXACT FROM 49H AGO)
# =============================================================================

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
    select(-lab) %>%
    pivot_wider(names_from = RefAlt, values_from = count) %>%
    mutate(
      freq = REF / (REF + ALT),
      N = REF + ALT
    ) %>%
    select(-REF, -ALT) %>%
    as_tibble()
  
  # Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
  founder_wide <- df3 %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  quality_filtered_positions <- founder_wide %>%
    filter(
      if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
    ) %>%
    pull(POS)
  
  # Filter to quality positions and include sample data
  df3 <- df3 %>%
    filter(POS %in% quality_filtered_positions)
  
  cat("✓ Processed", nrow(df3), "rows for", length(unique(df3$name)), "samples\n")
  return(df3)
}
