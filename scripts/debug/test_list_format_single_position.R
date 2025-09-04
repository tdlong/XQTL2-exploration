#!/usr/bin/env Rscript

# Test the actual list format haplotype estimation code on a single position
# This tests the estimate_haplotypes_list_format function from run_haplotype_estimation_list_format.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/debug/test_list_format_single_position.R <chr> <param_file> <output_dir> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
test_position <- as.numeric(args[4])

cat("=== TESTING LIST FORMAT HAPLOTYPE ESTIMATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Test position:", test_position, "\n\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("H cutoff:", h_cutoff, "\n\n")

# Define euchromatin boundaries (same as the main script)
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load observed data from REFALT files (same as the main script)
cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data file not found: ", refalt_file)
}

# Load REFALT data (same as the main script)
refalt_data <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(refalt_data), "rows\n")

# Transform to frequencies (same as the main script)
cat("Converting counts to frequencies...\n")
refalt_processed <- refalt_data %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    total_count = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

# Filter for high-quality SNPs (same as the main script)
cat("Filtering for high-quality SNPs...\n")
good_snps <- refalt_processed %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(total_count == 0),
    not_fixed = sum(total_count != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Get valid SNPs for euchromatin (same as the main script)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

# Filter to euchromatin and only the samples we actually processed (same as the main script)
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Source the actual functions from the main script
cat("Loading the actual list format functions...\n")

# Source the existing working functions first
source("scripts/production/haplotype_estimation_functions.R")

# Define the list format functions directly (extracted from the main script)
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                                           method = c("fixed", "adaptive"),
                                           window_size_bp = NULL,
                                           chr = "chr2R",
                                           verbose = 0) {
  
  method <- match.arg(method)
  
  if (verbose >= 1) {
    cat(sprintf("Processing pos: %s, sample: %s, method: %s\n", 
                format(pos, big.mark=","), sample_name, method))
  }
  
  # Validate inputs
  if (method == "fixed" && is.null(window_size_bp)) {
    stop("window_size_bp required for fixed method")
  }
  
  if (method == "fixed") {
    # FIXED WINDOW METHOD - use existing working code
    result <- estimate_haplotypes(pos, sample_name, df3, founders, h_cutoff,
                                "fixed", window_size_bp, chr, verbose)
    
    # Convert to list format
    return(create_list_result_from_existing(result, founders, NULL, NULL))
    
  } else {
    # ADAPTIVE WINDOW METHOD - modify existing code to capture groups and error matrix
    window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
    
    final_result <- NULL
    final_n_groups <- 0
    final_window_size <- window_sizes[1]
    final_wide_data <- NULL
    previous_n_groups <- 0
    final_groups <- NULL
    final_error_matrix <- NULL
    
    # CONSTRAINT ACCUMULATION (same as existing working code)
    accumulated_constraints <- NULL
    accumulated_constraint_values <- NULL
    
    for (window_idx in seq_along(window_sizes)) {
      window_size <- window_sizes[window_idx]
      window_start <- max(1, pos - window_size/2)
      window_end <- pos + window_size/2
      
      # Get SNPs in window (same as existing working code)
      window_data <- df3 %>%
        filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
      
      if (nrow(window_data) == 0) next
      
      wide_data <- window_data %>%
        select(POS, name, freq) %>%
        pivot_wider(names_from = name, values_from = freq)
      
      if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next
      
      # Get founder matrix and sample frequencies (same as existing working code)
      founder_matrix <- wide_data %>%
        select(all_of(founders)) %>%
        as.matrix()
      sample_freqs <- wide_data %>%
        pull(!!sample_name)
      
      complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
      founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
      sample_freqs_clean <- sample_freqs[complete_rows]
      
      if (nrow(founder_matrix_clean) < 10) next
      
      # Always update final window info (same as existing working code)
      final_window_size <- window_size
      final_wide_data <- wide_data
      
      # Hierarchical clustering (same as existing working code)
      founder_dist <- dist(t(founder_matrix_clean))
      hclust_result <- hclust(founder_dist, method = "complete")
      groups <- cutree(hclust_result, h = h_cutoff)
      n_groups <- length(unique(groups))
      
      # Check if clustering improved (same as existing working code)
      if (window_idx > 1 && n_groups <= previous_n_groups) {
        next  # No improvement, try larger window
      }
      
      previous_n_groups <- n_groups
      
      # Build constraint matrix (same as existing working code)
      n_founders <- ncol(founder_matrix_clean)
      E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      
      # Add accumulated constraints from previous windows (same as existing working code)
      if (!is.null(accumulated_constraints)) {
        E <- rbind(E, accumulated_constraints)
        F <- c(F, accumulated_constraint_values)
      }
      
      # Run LSEI with constraints AND fulloutput=TRUE (ONLY CHANGE!)
      tryCatch({
        result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                E = E, F = F, 
                                G = diag(n_founders), H = matrix(rep(0.0003, n_founders)),
                                fulloutput = TRUE)  # THIS IS THE ONLY CHANGE!
        
        if (result$IsError == 0) {
          # LSEI successful - capture the results (same as existing working code)
          final_result <- result
          final_n_groups <- n_groups
          final_groups <- groups
          names(final_groups) <- founders
          
          # Capture the error matrix (NEW!)
          if (!is.null(result$cov)) {
            final_error_matrix <- result$cov
            rownames(final_error_matrix) <- founders
            colnames(final_error_matrix) <- founders
          }
          
          # Accumulate constraints for next window (same as existing working code)
          current_constraints <- NULL
          current_constraint_values <- NULL
          
          unique_clusters <- unique(groups)
          for (cluster_id in unique_clusters) {
            cluster_founders <- which(groups == cluster_id)
            if (length(cluster_founders) > 1) {
              # Create constraint row for this group
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              
              # Calculate the actual group frequency from lsei result
              group_freq <- sum(result$X[cluster_founders])
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, group_freq)
            } else {
              # Single founder: lock their exact frequency
              founder_freq <- result$X[cluster_founders]
              
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, founder_freq)
            }
          }
          
          # Update accumulated constraints for next window (same as existing working code)
          if (!is.null(current_constraints)) {
            accumulated_constraints <- current_constraints
            accumulated_constraint_values <- current_constraint_values
          }
        }
      }, error = function(e) {
        # LSEI failed, continue to next window (same as existing working code)
      })
    }
    
    # Apply the correct rules for output (same as existing working code)
    if (!is.null(final_result)) {
      # LSEI was successful - get the results
      haplotype_freqs <- final_result$X
      names(haplotype_freqs) <- founders
      
      # Use the captured groups
      groups <- final_groups
      
      # Use the captured error matrix
      error_matrix <- final_error_matrix
      
      # Check if founders are distinguishable to set trust level (same as existing working code)
      if (final_n_groups == length(founders)) {
        estimate_OK <- 1  # Founders distinguishable
      } else {
        estimate_OK <- 0  # Founders NOT distinguishable
      }
      
      final_window_size <- final_window_size
      n_snps <- nrow(final_wide_data)
      
    } else {
      # Either insufficient SNPs OR LSEI failed/didn't converge (same as existing working code)
      haplotype_freqs <- rep(NA, length(founders))
      names(haplotype_freqs) <- founders
      groups <- rep(1, length(founders))
      names(groups) <- founders
      error_matrix <- matrix(NA, length(founders), length(founders))
      rownames(error_matrix) <- founders
      colnames(error_matrix) <- founders
      estimate_OK <- NA
      final_window_size <- window_sizes[1]
      n_snps <- 0
    }
    
    return(create_list_result(chr, pos, sample_name, method, final_window_size, n_snps, 
                            estimate_OK, haplotype_freqs, groups, error_matrix, founders))
  }
}

# Helper function to convert existing result to list format
create_list_result_from_existing <- function(existing_result, founders, groups, error_matrix) {
  # For fixed method, we don't have groups or error matrix from existing code
  if (is.null(groups)) {
    groups <- rep(1, length(founders))  # Default fallback
    names(groups) <- founders
  }
  if (is.null(error_matrix)) {
    error_matrix <- matrix(NA, length(founders), length(founders))
    rownames(error_matrix) <- founders
    colnames(error_matrix) <- founders
  }
  
  return(create_list_result(existing_result$chr, existing_result$pos, existing_result$sample, 
                          existing_result$method, existing_result$final_window_size, 
                          existing_result$n_snps, existing_result$estimate_OK, 
                          existing_result$haplotype_freqs, groups, error_matrix, founders))
}

# Create list format result
create_list_result <- function(chr, pos, sample_name, method, final_window_size, n_snps, 
                             estimate_OK, haplotype_freqs, groups, error_matrix, founders) {
  
  # Create the list format structure
  result <- list(
    CHROM = chr,
    pos = pos,
    sample = list(sample_name),
    Groups = list(list(groups)),
    Haps = list(list(haplotype_freqs)),
    Err = list(list(error_matrix)),
    Names = list(list(founders))
  )
  
  return(result)
}

# Test on first sample
test_sample <- names_in_bam[1]
cat("Testing on sample:", test_sample, "\n\n")

# Test the actual estimate_haplotypes_list_format function
cat("=== TESTING estimate_haplotypes_list_format FUNCTION ===\n")
result <- estimate_haplotypes_list_format(test_position, test_sample, observed_euchromatic, 
                                        founders, h_cutoff, "adaptive", 
                                        NULL, chr, verbose = 1)

cat("\n=== RESULTS ===\n")
cat("Result structure:\n")
str(result, max.level = 2)

cat("\nDetailed results:\n")
cat("CHROM:", result$CHROM, "\n")
cat("Position:", result$pos, "\n")
cat("Sample:", result$sample[[1]], "\n")

# Check groups
groups <- result$Groups[[1]][[1]]
cat("Groups:", paste(groups, collapse = ", "), "\n")
cat("Number of unique groups:", length(unique(groups)), "\n")
cat("Is 1:8?", all(sort(unique(groups)) == 1:8), "\n")

# Check haplotypes
haps <- result$Haps[[1]][[1]]
cat("Haplotypes:\n")
for (i in 1:length(haps)) {
  cat("  ", founders[i], ":", round(haps[i], 4), "\n")
}
cat("Sum:", round(sum(haps, na.rm = TRUE), 4), "\n")
cat("All NA?", all(is.na(haps)), "\n")

# Check error matrix
err <- result$Err[[1]][[1]]
cat("Error matrix dimensions:", nrow(err), "x", ncol(err), "\n")
cat("Error matrix all NA?", all(is.na(err)), "\n")
if (!all(is.na(err))) {
  cat("Error matrix has values:", sum(!is.na(err)), "/", length(err), "\n")
}

cat("\n=== TEST COMPLETE ===\n")
