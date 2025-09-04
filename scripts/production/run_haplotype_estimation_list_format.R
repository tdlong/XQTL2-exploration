#!/usr/bin/env Rscript

# Run Haplotype Estimation with List Format Output
# Focuses on adaptive_h4 and smooth_h4 methods
# Saves results to new hap_list_results directory

library(tidyverse)
library(limSolve)

# Source the new estimator functions
source("scripts/production/haplotype_estimation_functions.R")

# NEW HAPLOTYPE ESTIMATOR WITH LIST FORMAT OUTPUT
# This is a new function that captures groups and error matrices from the existing algorithm

estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                                           method = c("fixed", "adaptive"),
                                           window_size_bp = NULL,
                                           chr = "chr2R",
                                           verbose = 0) {
  
  method <- match.arg(method)
  
  # Initialize result variables
  estimate_OK <- NA
  haplotype_freqs <- rep(NA, length(founders))
  names(haplotype_freqs) <- founders
  error_matrix <- matrix(NA, length(founders), length(founders))
  rownames(error_matrix) <- founders
  colnames(error_matrix) <- founders
  groups <- rep(NA, length(founders))
  names(groups) <- founders
  final_window_size <- NA
  n_snps <- 0
  
  if (method == "fixed") {
    # FIXED WINDOW METHOD
    window_start <- max(1, pos - window_size_bp/2)
    window_end <- pos + window_size_bp/2
    
    # Get SNPs in window
    window_data <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    if (nrow(window_data) == 0) {
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                              rep(NA, length(founders)), rep(1, length(founders)), 
                              matrix(NA, length(founders), length(founders)), founders))
    }
    
    wide_data <- window_data %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) {
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                              rep(NA, length(founders)), rep(1, length(founders)), 
                              matrix(NA, length(founders), length(founders)), founders))
    }
    
    # Get founder matrix and sample frequencies
    founder_matrix_clean <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    sample_freqs_clean <- wide_data %>%
      pull(!!sample_name)
    
    complete_rows <- complete.cases(founder_matrix_clean) & !is.na(sample_freqs_clean)
    founder_matrix_clean <- founder_matrix_clean[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs_clean[complete_rows]
    
    if (nrow(founder_matrix_clean) < 10) {
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                              rep(NA, length(founders)), rep(1, length(founders)), 
                              matrix(NA, length(founders), length(founders)), founders))
    }
    
    final_window_size <- window_size_bp
    n_snps <- nrow(wide_data)
    
    # Run LSEI with error matrix capture
    tryCatch({
      E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      G <- diag(length(founders))  # Non-negativity constraints
      H <- matrix(rep(0.0003, length(founders)))  # Lower bound
      
      # Call lsei with fulloutput=TRUE to get error matrix
      lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                   E = E, F = F, G = G, H = H, fulloutput = TRUE)
      
      if (lsei_result$IsError == 0) {
        # LSEI successful - get frequencies
        haplotype_freqs <- lsei_result$X
        names(haplotype_freqs) <- founders
        
        # Capture the error matrix
        if (!is.null(lsei_result$cov)) {
          error_matrix <- lsei_result$cov
          rownames(error_matrix) <- founders
          colnames(error_matrix) <- founders
        }
        
        # Get clustering groups (this is the key part that was missing!)
        distances <- dist(t(founder_matrix_clean))
        hclust_result <- hclust(distances, method = "complete")
        groups <- cutree(hclust_result, h = h_cutoff)
        names(groups) <- founders
        
        # Determine estimate_OK based on distinguishability
        n_groups <- length(unique(groups))
        estimate_OK <- ifelse(n_groups == length(founders), 1, 0)
        
      } else {
        # LSEI failed
        estimate_OK <- 0
        groups <- rep(1, length(founders))
        names(groups) <- founders
      }
    }, error = function(e) {
      # Catastrophic LSEI error
      estimate_OK <- 0
      groups <- rep(1, length(founders))
      names(groups) <- founders
    })
    
  } else {
    # ADAPTIVE WINDOW METHOD - implement the full adaptive algorithm with groups and error matrix capture
    window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
    
    final_result <- NULL
    final_n_groups <- 0
    final_window_size <- window_sizes[1]
    final_wide_data <- NULL
    previous_n_groups <- 0
    final_groups <- NULL
    final_error_matrix <- NULL
    
    # CONSTRAINT ACCUMULATION (core of adaptive algorithm)
    accumulated_constraints <- NULL
    accumulated_constraint_values <- NULL
    
    for (window_idx in seq_along(window_sizes)) {
      window_size <- window_sizes[window_idx]
      window_start <- max(1, pos - window_size/2)
      window_end <- pos + window_size/2
      
      # Get SNPs in window
      window_data <- df3 %>%
        filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
      
      if (nrow(window_data) == 0) next
      
      wide_data <- window_data %>%
        select(POS, name, freq) %>%
        pivot_wider(names_from = name, values_from = freq)
      
      if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next
      
      # Get founder matrix and sample frequencies
      founder_matrix <- wide_data %>%
        select(all_of(founders)) %>%
        as.matrix()
      sample_freqs <- wide_data %>%
        pull(!!sample_name)
      
      complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
      founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
      sample_freqs_clean <- sample_freqs[complete_rows]
      
      if (nrow(founder_matrix_clean) < 10) next
      
      # Always update final window info (this is the last window we actually tried)
      final_window_size <- window_size
      final_wide_data <- wide_data
      
      # Hierarchical clustering
      founder_dist <- dist(t(founder_matrix_clean))
      hclust_result <- hclust(founder_dist, method = "complete")
      groups <- cutree(hclust_result, h = h_cutoff)
      n_groups <- length(unique(groups))
      
      # Check if clustering improved
      if (window_idx > 1 && n_groups <= previous_n_groups) {
        next  # No improvement, try larger window
      }
      
      previous_n_groups <- n_groups
      
      # Build constraint matrix with accumulated constraints
      n_founders <- ncol(founder_matrix_clean)
      E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      
      # Add accumulated constraints from previous windows
      if (!is.null(accumulated_constraints)) {
        E <- rbind(E, accumulated_constraints)
        F <- c(F, accumulated_constraint_values)
      }
      
      # Run LSEI with constraints AND fulloutput=TRUE to get error matrix
      tryCatch({
        result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                E = E, F = F, 
                                G = diag(n_founders), H = matrix(rep(0.0003, n_founders)),
                                fulloutput = TRUE)  # THIS IS THE KEY - get error matrix
        
        if (result$IsError == 0) {
          # LSEI successful - capture the results
          final_result <- result
          final_n_groups <- n_groups
          final_groups <- groups
          names(final_groups) <- founders
          
          # Capture the error matrix
          if (!is.null(result$cov)) {
            final_error_matrix <- result$cov
            rownames(final_error_matrix) <- founders
            colnames(final_error_matrix) <- founders
          }
          
          # Accumulate constraints for next window (same as original algorithm)
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
          
          # Update accumulated constraints for next window
          if (!is.null(current_constraints)) {
            accumulated_constraints <- current_constraints
            accumulated_constraint_values <- current_constraint_values
          }
        }
      }, error = function(e) {
        # LSEI failed, continue to next window
      })
    }
    
    # Apply the correct rules for output
    if (!is.null(final_result)) {
      # LSEI was successful - get the results
      haplotype_freqs <- final_result$X
      names(haplotype_freqs) <- founders
      
      # Use the captured groups
      groups <- final_groups
      
      # Use the captured error matrix
      error_matrix <- final_error_matrix
      
      # The groups vector IS the information - no need for estimate_OK
      # Groups will be 1:8 if all founders distinguishable, or something else if not
      estimate_OK <- NA  # We don't use estimate_OK anymore - groups contains the info
      
      final_window_size <- final_window_size
      n_snps <- nrow(final_wide_data)
      
    } else {
      # Either insufficient SNPs OR LSEI failed/didn't converge
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
  }
  
  return(create_list_result(chr, pos, sample_name, method, final_window_size, n_snps, 
                          estimate_OK, haplotype_freqs, groups, error_matrix, founders))
}

# Helper function to create the list format result
create_list_result <- function(chr, pos, sample_name, method, window_size, n_snps, 
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

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/production/run_haplotype_estimation_list_format.R <chr> <param_file> <output_dir> <method>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
method <- args[4]

# Validate method
if (!method %in% c("adaptive_h4", "smooth_h4")) {
  stop("Method must be 'adaptive_h4' or 'smooth_h4'")
}

cat("=== HAPLOTYPE ESTIMATION WITH LIST FORMAT ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Method:", method, "\n\n")

# Load parameters
source(param_file)

# Define euchromatin boundaries (same as existing pipeline)
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

# Get boundaries for this chromosome
if (!chr %in% names(euchromatin_boundaries)) {
  stop("Invalid chromosome: ", chr, ". Valid chromosomes: ", paste(names(euchromatin_boundaries), collapse = ", "))
}

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region for", chr, ":", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load observed data from REFALT files (same as existing pipeline)
cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data file not found: ", refalt_file)
}

# Load REFALT data
refalt_data <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(refalt_data), "rows\n")

# Transform to frequencies (same as existing pipeline)
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

# Filter for high-quality SNPs (same as existing pipeline)
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

# Get valid SNPs for evaluation (euchromatin only)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

# Filter to euchromatin and only the samples we actually processed
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)  # Only process samples defined in parameter file

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Create new results directory
list_results_dir <- file.path(output_dir, "list_results")
dir.create(list_results_dir, showWarnings = FALSE, recursive = TRUE)

# Define positions to process
positions <- seq(euchromatin_start, euchromatin_end, by = 10000)
samples <- names_in_bam

cat("Processing", length(positions), "positions for", length(samples), "samples\n")
cat("Method:", method, "\n\n")

if (method == "adaptive_h4") {
  # Run adaptive_h4 with list format
  cat("Running adaptive_h4 estimation...\n")
  
  list_format_data <- positions %>%
    map_dfr(function(pos) {
      if (pos %% 100000 == 0) {
        cat("Processing position:", format(pos, big.mark=","), "\n")
      }
      
      # Process each sample at this position
      sample_results <- samples %>%
        map(function(sample_name) {
          result <- estimate_haplotypes_list_format(pos, sample_name, observed_euchromatic, 
                                                  founders, h_cutoff, "adaptive", 
                                                  NULL, chr, verbose = 0)
          return(result)
        }) %>%
        compact()  # Remove NULL results
      
      if (length(sample_results) == 0) {
        return(NULL)
      }
      
      # Combine results for all samples at this position
      all_samples <- unlist(map(sample_results, "sample"))
      all_groups <- map(sample_results, ~ .x$Groups[[1]][[1]])
      all_haps <- map(sample_results, ~ .x$Haps[[1]][[1]])
      all_err <- map(sample_results, ~ .x$Err[[1]][[1]])
      
      # Create the combined list format structure for this position
      tibble(
        CHROM = chr,
        pos = pos,
        sample = list(all_samples),
        Groups = list(all_groups),
        Haps = list(all_haps),
        Err = list(all_err),
        Names = list(rep(list(founders), length(all_samples)))
      )
    }) %>%
    filter(!is.null(pos))  # Remove NULL results
  
  # Save adaptive_h4 results
  output_file <- file.path(list_results_dir, paste0("adaptive_h4_list_format_", chr, ".RDS"))
  saveRDS(list_format_data, output_file)
  cat("✓ Adaptive_h4 list format results saved:", nrow(list_format_data), "rows\n")
  cat("Output file:", output_file, "\n")
  
} else if (method == "smooth_h4") {
  # Run smooth_h4 with list format
  cat("Running smooth_h4 estimation...\n")
  
  # First, load the adaptive_h4 results to smooth over
  adaptive_file <- file.path(list_results_dir, paste0("adaptive_h4_list_format_", chr, ".RDS"))
  if (!file.exists(adaptive_file)) {
    stop("Adaptive_h4 results not found. Run adaptive_h4 first: ", adaptive_file)
  }
  
  adaptive_data <- readRDS(adaptive_file)
  cat("✓ Loaded adaptive_h4 results:", nrow(adaptive_data), "positions\n")
  
  # Create smooth_h4 by averaging adaptive_h4 results
  # For smooth_h4, we use a 21-position sliding window (matching original implementation)
  window_size <- 21
  
  # First, we need to determine the quality of each position based on groups
  # We'll check groups directly in the smooth_h4 logic - no need for position_ok
  adaptive_data_with_quality <- adaptive_data
  
  # Apply 21-position sliding window smoothing
  smooth_data <- adaptive_data_with_quality %>%
    arrange(pos) %>%
    mutate(
                 # Calculate quality count: how many positions in 21-position window have 8 groups?
           quality_count = map_dbl(seq_len(n()), function(i) {
             start_idx <- max(1, i - 10)  # 10 positions before
             end_idx <- min(n(), i + 10)  # 10 positions after
             window_positions <- start_idx:end_idx
             
             # Check if each position has 8 distinguishable groups (1:8)
             positions_with_8_groups <- map_lgl(window_positions, function(j) {
               groups_at_pos <- adaptive_data_with_quality$Groups[[j]][[1]]
               if (is.null(groups_at_pos) || any(is.na(groups_at_pos))) return(FALSE)
               # Check if all 8 founders are in different groups (1:8)
               length(unique(groups_at_pos)) == 8 && all(sort(unique(groups_at_pos)) == 1:8)
             })
             
             sum(positions_with_8_groups, na.rm = TRUE)
           }),
      
      # New estimate_OK: OK if at least 17 out of 21 positions are OK
      new_estimate_ok = quality_count >= 17,
      
      # For smooth_h4, groups depend on the quality
      Groups = list(map(seq_along(sample[[1]]), function(i) {
        if (new_estimate_ok) {
          # If smooth estimate is OK, use 1:8 (all founders distinguishable)
          return(1:8)
        } else {
          # If smooth estimate is not OK, use all 1s (fallback)
          return(rep(1, length(founders)))
        }
      })),
      
      # Average the haplotype frequencies using 21-position sliding window
      Haps = list(map(seq_along(sample[[1]]), function(i) {
        pos_idx <- which(adaptive_data_with_quality$pos == pos)
        start_idx <- max(1, pos_idx - 10)  # 10 positions before
        end_idx <- min(nrow(adaptive_data_with_quality), pos_idx + 10)  # 10 positions after
        
                     # Get haplotype estimates for this sample across the 21-position window
             window_haps <- map(start_idx:end_idx, function(j) {
               adaptive_data_with_quality$Haps[[j]][[1]][[i]]
             })
             
             # Only use positions where the original estimate had 8 groups
             valid_positions <- map_lgl(start_idx:end_idx, function(j) {
               groups_at_pos <- adaptive_data_with_quality$Groups[[j]][[1]]
               if (is.null(groups_at_pos) || any(is.na(groups_at_pos))) return(FALSE)
               # Check if all 8 founders are in different groups (1:8)
               length(unique(groups_at_pos)) == 8 && all(sort(unique(groups_at_pos)) == 1:8)
             })
        
        # Filter to only OK positions (not all 21 positions)
        valid_haps <- window_haps[valid_positions & map_lgl(window_haps, ~ !any(is.na(.x)))]
        
        if (new_estimate_ok && length(valid_haps) > 0) {
          # Average over only the OK positions (not all 21 positions)
          avg_haps <- reduce(valid_haps, `+`) / length(valid_haps)
          # Normalize so they sum to 1
          avg_haps <- avg_haps / sum(avg_haps)
          names(avg_haps) <- founders
          return(avg_haps)
        } else {
          # If <17 positions are OK, return NAs
          return(set_names(rep(NA, length(founders)), founders))
        }
      })),
      
      # Average the error matrices using 21-position sliding window
      Err = list(map(seq_along(sample[[1]]), function(i) {
        pos_idx <- which(adaptive_data_with_quality$pos == pos)
        start_idx <- max(1, pos_idx - 10)  # 10 positions before
        end_idx <- min(nrow(adaptive_data_with_quality), pos_idx + 10)  # 10 positions after
        
                     # Get error matrices for this sample across the 21-position window
             window_errs <- map(start_idx:end_idx, function(j) {
               adaptive_data_with_quality$Err[[j]][[1]][[i]]
             })
             
             # Only use positions where the original estimate had 8 groups
             valid_positions <- map_lgl(start_idx:end_idx, function(j) {
               groups_at_pos <- adaptive_data_with_quality$Groups[[j]][[1]]
               if (is.null(groups_at_pos) || any(is.na(groups_at_pos))) return(FALSE)
               # Check if all 8 founders are in different groups (1:8)
               length(unique(groups_at_pos)) == 8 && all(sort(unique(groups_at_pos)) == 1:8)
             })
        
        # Filter to only OK positions (not all 21 positions)
        valid_errs <- window_errs[valid_positions & map_lgl(window_errs, ~ !any(is.na(.x)))]
        
        if (new_estimate_ok && length(valid_errs) > 0) {
          # Average over only the OK positions (not all 21 positions)
          avg_err <- reduce(valid_errs, `+`) / length(valid_errs)
          rownames(avg_err) <- founders
          colnames(avg_err) <- founders
          return(avg_err)
        } else {
          # If <17 positions are OK, return NAs
          return(matrix(NA, length(founders), length(founders), 
                       dimnames = list(founders, founders)))
        }
      }))
    ) %>%
             select(-quality_count, -new_estimate_ok)  # Clean up temporary columns
  
  # Save smooth_h4 results
  output_file <- file.path(list_results_dir, paste0("smooth_h4_list_format_", chr, ".RDS"))
  saveRDS(smooth_data, output_file)
  cat("✓ Smooth_h4 list format results saved:", nrow(smooth_data), "rows\n")
  cat("Output file:", output_file, "\n")
}

cat("\n=== HAPLOTYPE ESTIMATION COMPLETE ===\n")
cat("Results saved to:", list_results_dir, "\n")
