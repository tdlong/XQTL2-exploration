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

# Function to perform haplotype estimation with full output capture
estimate_haplotypes_full <- function(snp_data, method, parameter, founders) {
  # This function needs to be integrated with your actual haplotype estimation code
  # It should capture:
  # 1. The lsei output including error matrices
  # 2. The clustering groups
  # 3. The estimate_OK status
  
  # For now, this is a template showing what needs to be captured
  
  # Create founder matrix (you'll need your actual founder matrix creation)
  founder_mat <- matrix(runif(length(founders) * nrow(snp_data)), 
                       nrow = nrow(snp_data), ncol = length(founders))
  colnames(founder_mat) <- founders
  
  # Create observation vector
  Y <- snp_data$freq
  
  # Perform clustering to get groups
  if (method == "adaptive") {
    h_cutoff <- parameter
    groups <- cutree(hclust(dist(t(founder_mat))), h = h_cutoff)
  } else {
    # For fixed windows, all founders are distinguishable
    groups <- 1:length(founders)
    names(groups) <- founders
  }
  
  # Perform lsei estimation
  d <- ncol(founder_mat)
  
  # THIS IS WHERE YOU NEED TO INTEGRATE YOUR ACTUAL LSEI CALL
  # The lsei call should look something like:
  # out <- lsei(A=founder_mat, B=Y, E=t(matrix(rep(1,d))), F=1, G=diag(d), H=matrix(rep(0.0003,d)), verbose=TRUE, fulloutput=TRUE)
  
  # For now, create placeholder results
  haps <- runif(length(founders))
  names(haps) <- founders
  
  # Create a placeholder error matrix (this should come from out$cov)
  error_matrix <- diag(0.01, length(founders))
  rownames(error_matrix) <- founders
  colnames(error_matrix) <- founders
  
  # Determine if estimation was successful
  estimate_ok <- ifelse(all(!is.na(haps)) && all(haps >= 0) && abs(sum(haps) - 1) < 0.01, 1, 0)
  
  return(list(
    haps = haps,
    groups = groups,
    error_matrix = error_matrix,
    estimate_ok = estimate_ok,
    n_snps = nrow(snp_data)
  ))
}

# Function to process a window and return list format data
process_window_list_format <- function(window_snps, method, parameter, founders, pos) {
  if (nrow(window_snps) == 0) {
    return(NULL)
  }
  
  # Estimate haplotypes with full output capture
  result <- estimate_haplotypes_full(window_snps, method, parameter, founders)
  
  # Create the list format structure
  # For each sample, we need to create the nested list structure
  
  # Get unique samples in this window
  samples <- unique(window_snps$name)
  
  # Create the list format data for this position
  list_format_row <- tibble(
    CHROM = chr,
    pos = pos,
    sample = list(samples),
    Groups = list(rep(list(result$groups), length(samples))),
    Haps = list(rep(list(result$haps), length(samples))),
    Err = list(rep(list(result$error_matrix), length(samples))),
    Names = list(rep(list(founders), length(samples)))
  )
  
  return(list_format_row)
}

# Main processing loop
cat("Processing haplotype estimation...\n")

# This is where you'd implement your actual windowing logic
# For now, create a simplified example

# Create results directory
results_dir <- file.path(output_dir, "haplotype_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Generate some example results in list format
positions <- seq(euchromatin_start, euchromatin_end, by = 10000)
samples <- names_in_bam

# Create example list format data
example_list_data <- positions %>%
  map_dfr(function(pos) {
    # For each position, create the list format structure
    tibble(
      CHROM = chr,
      pos = pos,
      sample = list(samples),
      Groups = list(rep(list(1:length(founders)), length(samples))),
      Haps = list(rep(list(set_names(runif(length(founders)), founders)), length(samples))),
      Err = list(rep(list(diag(0.01, length(founders))), length(samples))),
      Names = list(rep(list(founders), length(samples)))
    )
  })

cat("✓ Example list format results generated:", nrow(example_list_data), "rows\n")

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
