Great.  #!/usr/bin/env Rscript

# Test both the working function and the new list format function on the same position
# This ensures the new function gives the same results as the working function

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

cat("=== TESTING BOTH WORKING AND NEW LIST FORMAT FUNCTIONS ===\n")
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

# Load REFALT data (same as existing working code)
refalt_data <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(refalt_data), "rows\n")

# Transform to frequencies (same as existing working code)
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

# Filter for high-quality SNPs (same as existing working code)
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

# Get valid SNPs for euchromatin (same as existing working code)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

# Filter to euchromatin and only the samples we actually processed (same as existing working code)
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Source the existing working functions
source("scripts/production/haplotype_estimation_functions.R")

# Define the new function (copied from run_haplotype_estimation_list_format.R)
estimate_haplotypes_with_groups <- function(pos, sample_name, df3, founders, h_cutoff, chr, verbose = 0) {
  
  if (verbose >= 1) {
    cat(sprintf("Processing pos: %s, sample: %s, method: adaptive\n", 
                format(pos, big.mark=","), sample_name))
  }
  
  # ADAPTIVE WINDOW METHOD - copied from working code
  window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
  
  final_result <- NULL
  final_n_groups <- 0
  final_window_size <- window_sizes[1]
  final_wide_data <- NULL
  previous_n_groups <- 0
  final_groups <- NULL
  final_error_matrix <- NULL
  
  # CONSTRAINT ACCUMULATION (copied from working code)
  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL
  
  for (window_idx in seq_along(window_sizes)) {
    window_size <- window_sizes[window_idx]
    window_start <- max(1, pos - window_size/2)
    window_end <- pos + window_size/2
    
    # Get SNPs in window (copied from working code)
    window_data <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    if (nrow(window_data) == 0) next
    
    wide_data <- window_data %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next
    
    # Get founder matrix and sample frequencies (copied from working code)
    founder_matrix <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    sample_freqs <- wide_data %>%
      pull(!!sample_name)
    
    complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    
    if (nrow(founder_matrix_clean) < 10) next
    
    # Always update final window info (copied from working code)
    final_window_size <- window_size
    final_wide_data <- wide_data
    
    # Hierarchical clustering (copied from working code)
    founder_dist <- dist(t(founder_matrix_clean))
    hclust_result <- hclust(founder_dist, method = "complete")
    groups <- cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    
    # Check if clustering improved (copied from working code)
    if (window_idx > 1 && n_groups <= previous_n_groups) {
      next  # No improvement, try larger window
    }
    
    previous_n_groups <- n_groups
    
    # Build constraint matrix (copied from working code)
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    
    # Add accumulated constraints from previous windows (copied from working code)
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
        # LSEI successful - capture the results (copied from working code)
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
        
        # Accumulate constraints for next window (copied from working code)
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
        
        # Update accumulated constraints for next window (copied from working code)
        if (!is.null(current_constraints)) {
          accumulated_constraints <- current_constraints
          accumulated_constraint_values <- current_constraint_values
        }
      }
    }, error = function(e) {
      # LSEI failed, continue to next window (copied from working code)
    })
  }
  
  # Apply the correct rules for output (copied from working code)
  if (!is.null(final_result)) {
    # LSEI was successful - get the results
    haplotype_freqs <- final_result$X
    names(haplotype_freqs) <- founders
    
    # Use the captured groups
    groups <- final_groups
    
    # Use the captured error matrix
    error_matrix <- final_error_matrix
    
    # Check if founders are distinguishable to set trust level (copied from working code)
    if (final_n_groups == length(founders)) {
      estimate_OK <- 1  # Founders distinguishable
    } else {
      estimate_OK <- 0  # Founders NOT distinguishable
    }
    
    final_window_size <- final_window_size
    n_snps <- nrow(final_wide_data)
    
  } else {
    # Either insufficient SNPs OR LSEI failed/didn't converge (copied from working code)
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
  
  # Return in the same format as the working function, plus groups and error matrix
  return(list(
    chr = chr,
    pos = pos,
    sample = sample_name,
    method = "adaptive",
    final_window_size = final_window_size,
    n_snps = n_snps,
    estimate_OK = estimate_OK,
    haplotype_freqs = haplotype_freqs,
    groups = groups,
    error_matrix = error_matrix
  ))
}

# Test both functions on the same position and sample
test_sample <- names_in_bam[1]  # Use first sample
cat("=== TESTING BOTH FUNCTIONS ON SAME POSITION ===\n")
cat("Test position:", test_position, "\n")
cat("Test sample:", test_sample, "\n\n")

# Test the existing working function first
cat("=== TESTING EXISTING WORKING FUNCTION ===\n")
working_result <- estimate_haplotypes(test_position, test_sample, observed_euchromatic, 
                                    founders, h_cutoff, "adaptive", 
                                    NULL, chr, verbose = 2)

cat("\nWorking function results:\n")
cat("Estimate OK:", working_result$estimate_OK, "\n")
cat("Final window size:", working_result$final_window_size, "\n")
cat("Number of SNPs:", working_result$n_snps, "\n")
cat("Haplotype frequencies:\n")
print(working_result$haplotype_freqs)

# Test the new function with groups and error matrix
cat("\n=== TESTING NEW FUNCTION WITH GROUPS AND ERROR MATRIX ===\n")
new_result <- estimate_haplotypes_with_groups(test_position, test_sample, observed_euchromatic, 
                                             founders, h_cutoff, chr, verbose = 2)

cat("\nNew function results:\n")
cat("Estimate OK:", new_result$estimate_OK, "\n")
cat("Final window size:", new_result$final_window_size, "\n")
cat("Number of SNPs:", new_result$n_snps, "\n")
cat("Haplotype frequencies:\n")
print(new_result$haplotype_freqs)
cat("Groups:\n")
print(new_result$groups)
cat("Error matrix dimensions:", dim(new_result$error_matrix), "\n")
if (!all(is.na(new_result$error_matrix))) {
  cat("Error matrix (first few values):\n")
  print(new_result$error_matrix[1:3, 1:3])
} else {
  cat("Error matrix: All NAs\n")
}

# Compare the results
cat("\n=== COMPARISON ===\n")
cat("Estimate OK match:", working_result$estimate_OK == new_result$estimate_OK, "\n")
cat("Window size match:", working_result$final_window_size == new_result$final_window_size, "\n")
cat("SNP count match:", working_result$n_snps == new_result$n_snps, "\n")

# Compare haplotype frequencies
if (!any(is.na(working_result$haplotype_freqs)) && !any(is.na(new_result$haplotype_freqs))) {
  freq_diff <- abs(working_result$haplotype_freqs - new_result$haplotype_freqs)
  cat("Max haplotype frequency difference:", max(freq_diff), "\n")
  cat("Haplotype frequencies match (within 1e-10):", all(freq_diff < 1e-10), "\n")
} else {
  cat("Haplotype frequencies: One or both contain NAs\n")
}

cat("\n=== TEST COMPLETE ===\n")