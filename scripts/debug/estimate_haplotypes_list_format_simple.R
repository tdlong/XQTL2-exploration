#!/usr/bin/env Rscript

# Simple wrapper that calls the working function and converts to list format
# This just changes the return format, doesn't rewrite the working logic

source("scripts/production/haplotype_estimation_functions.R")

#' Estimate haplotypes using list format (wrapper around working function)
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff, method, window_size_bp, chr, verbose = 0) {
  
  # Call the working function exactly as it is
  result <- estimate_haplotypes(
    pos = pos,
    sample_name = sample_name, 
    df3 = df3,
    founders = founders,
    h_cutoff = h_cutoff,
    method = method,
    window_size_bp = window_size_bp,
    chr = chr,
    verbose = verbose
  )
  
  # Convert the working result to list format
  if (result$estimate_OK) {
    # Extract haplotype frequencies for this sample
    haplotype_freqs <- numeric(length(founders))
    for (i in seq_along(founders)) {
      haplotype_freqs[i] <- result[[founders[i]]]
    }
    
    # For now, create placeholder groups and error matrix
    # TODO: Extract these from the working function's internal results
    groups <- rep(1, length(founders))  # Placeholder
    error_matrix <- matrix(NA, length(founders), length(founders))  # Placeholder
    
    return(list(
      CHROM = chr,
      pos = pos,
      sample = sample_name,
      Groups = list(groups),
      Haps = list(haplotype_freqs),
      Err = list(error_matrix),
      Names = list(founders)
    ))
  } else {
    # Return empty result for failed estimation
    return(list(
      CHROM = chr,
      pos = pos,
      sample = sample_name,
      Groups = list(rep(1, length(founders))),
      Haps = list(rep(NA, length(founders))),
      Err = list(matrix(NA, length(founders), length(founders))),
      Names = list(founders)
    ))
  }
}


