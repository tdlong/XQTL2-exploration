#!/usr/bin/env Rscript

# New haplotype estimation function that outputs list format directly
# This function captures all required information:
# 1. Actual error matrices from lsei output
# 2. Meaningful groups from clustering
# 3. Proper handling of estimate_OK cases
# 4. Direct output in the desired list format

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript haplotype_estimation_list_format.R <chr> <method> <parameter> <output_dir> <param_file>")
}

chr <- args[1]
method <- args[2]
parameter <- as.numeric(args[3])
output_dir <- args[4]
param_file <- args[5]

cat("=== HAPLOTYPE ESTIMATION (LIST FORMAT) ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Load parameters
source(param_file)

# Define euchromatin boundaries
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]

cat("Euchromatin boundaries:", euchromatin_start, "to", euchromatin_end, "bp\n\n")

# Load REFALT data
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data not found: ", refalt_file)
}

cat("Loading REFALT data...\n")
refalt_data <- read.table(refalt_file, header = TRUE)

# Transform to frequencies
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

# Filter for high-quality SNPs
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

# Get valid SNPs for euchromatin
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

cat("✓ Valid euchromatic SNPs:", nrow(observed_data %>% distinct(CHROM, POS)), "\n")

# Filter to euchromatin and only the samples we actually processed
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

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
    # ADAPTIVE WINDOW METHOD - simplified for now
    # This would need the full adaptive logic from the existing code
    estimate_OK <- 0
    groups <- rep(1, length(founders))
    names(groups) <- founders
    final_window_size <- 100000  # placeholder
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

# The main processing loop uses the new estimate_haplotypes_list_format function directly

# Main processing loop
cat("Processing haplotype estimation...\n")

# Create results directory
results_dir <- file.path(output_dir, "haplotype_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Define positions to process
positions <- seq(euchromatin_start, euchromatin_end, by = 10000)
samples <- names_in_bam

# Process each position using the new estimator
list_format_data <- positions %>%
  map_dfr(function(pos) {
    # Process each sample at this position
    sample_results <- samples %>%
      map(function(sample_name) {
        result <- estimate_haplotypes_list_format(pos, sample_name, observed_euchromatic, 
                                                founders, parameter, method, 
                                                ifelse(method == "fixed", parameter * 1000, NULL),
                                                chr, verbose = 0)
        return(result)
      }) %>%
      compact()  # Remove NULL results
    
    if (length(sample_results) == 0) {
      return(NULL)
    }
    
    # Combine results for all samples at this position
    # Each sample result is already in the correct list format
    # We need to combine them into a single row per position
    
    # Extract the data from each sample result
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

cat("✓ List format results generated:", nrow(list_format_data), "rows\n")

# Save results
if (method == "fixed") {
  output_file <- file.path(results_dir, paste0("fixed_window_", parameter, "kb_list_format_", chr, ".RDS"))
} else if (method == "adaptive") {
  output_file <- file.path(results_dir, paste0("adaptive_window_h", parameter, "_list_format_", chr, ".RDS"))
} else {
  stop("Unknown method: ", method)
}

saveRDS(example_list_data, output_file)

cat("✓ Results saved to:", basename(output_file), "\n")
cat("File size:", round(file.size(output_file) / 1024^2, 2), "MB\n")

cat("\n=== HAPLOTYPE ESTIMATION COMPLETE ===\n")
cat("This is a template script that shows the structure needed.\n")
cat("You need to integrate this with your actual haplotype estimation code.\n")
cat("Key points:\n")
cat("1. Capture the actual lsei output including error matrices\n")
cat("2. Store the clustering groups from cutree\n")
cat("3. Output directly in the list format structure\n")
cat("4. Handle estimate_OK cases properly (groups=all 1s, haps=NA, err=NA)\n")
cat("5. For smooth_h4, average over the underlying adaptive_h4 results\n")
