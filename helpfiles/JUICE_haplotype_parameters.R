# JUICE Haplotype Parameters
# This file defines the parameters for haplotype estimation in the JUICE dataset

# Founder samples (must match REFALT founder columns)
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Step size for scanning positions along chromosome
step <- 10000  # 10kb steps

# Window size comes from command line; keep commented here
# size <- 50000

# NOTE: In adaptive mode, h_cutoff is provided via CLI/SLURM (parameter arg).
# Do not define h_cutoff here to avoid ambiguity.

# Samples to process (select subset from large REFALT files)
# Full 12-sample experiment
names_in_bam <- c("AJ_1_1", "AJ_2_1", "AJ_3_1", "AJ_4_1", "AJ_5_1", "AJ_6_1", 
                  "GJ_1_1", "GJ_2_1", "GJ_3_1", "GJ_4_1", "GJ_5_1", "GJ_6_1")
