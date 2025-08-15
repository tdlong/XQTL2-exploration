# JUICE Haplotype Parameters
# This file defines the parameters for haplotype estimation in the JUICE dataset

# Founder samples (these should match the column names in your REFALT data)
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Non-founder samples (these will be estimated)
# Note: These will be automatically detected from the REFALT data
# but you can specify them here if needed

# Quality thresholds
min_freq <- 0.03      # Minimum frequency for informative SNPs
max_freq <- 0.97      # Maximum frequency for informative SNPs
min_total_count <- 0  # Minimum total count for SNPs

# Window parameters
min_window_snps <- 10  # Minimum SNPs required in a window
max_window_size <- 500000  # Maximum window size in bp

# Clustering parameters
h_cutoff <- 6  # Default h_cutoff for hierarchical clustering

# Output parameters
output_dir <- "process/JUICE/haplotype_results"
