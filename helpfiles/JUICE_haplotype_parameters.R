# JUICE Haplotype Parameters
# This file defines the parameters for haplotype estimation in the JUICE dataset

# Founder samples (these should match the column names in your REFALT data)
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Step size for scanning positions along chromosome
step <- 10000  # 10kb steps (increased from 5kb for speed)

# Window size - COMMENTED OUT, comes from command line
# size <- 50000  # Window size comes from command line parameter

# Clustering parameters
h_cutoff <- 2.5  # Default h_cutoff for fixed window hierarchical clustering

# Samples to process (allows selecting subset from large REFALT files)
names_in_bam <- c("AJ_1_1", "AJ_2_1", "AJ_3_1", "GJ_1_1", "GJ_2_1", "GJ_3_1")
