# ZINC2 Haplotype Parameters
# This file defines the parameters for haplotype estimation in the ZINC2 dataset

# Founder samples (must match REFALT founder columns)
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Step size for scanning positions along chromosome
step <- 10000  # 10kb steps

# Window size comes from command line; keep commented here
# size <- 50000

# NOTE: In adaptive mode, h_cutoff is provided via CLI/SLURM (parameter arg).
# Do not define h_cutoff here to avoid ambiguity.

# Samples to process (select subset from large REFALT files)
# Full 60-sample experiment
names_in_bam <- c(
  "Rep01_W_F","Rep01_W_M","Rep01_Z_F","Rep01_Z_M",
  "Rep02_W_F","Rep02_W_M","Rep02_Z_F","Rep02_Z_M",
  "Rep03_W_F","Rep03_W_M","Rep03_Z_F","Rep03_Z_M",
  "Rep04_W_F","Rep04_W_M","Rep04_Z_F","Rep04_Z_M",
  "Rep05_W_F","Rep05_W_M","Rep05_Z_F","Rep05_Z_M",
  "Rep06_W_F","Rep06_W_M","Rep06_Z_F","Rep06_Z_M",
  "Rep07_W_F","Rep07_W_M","Rep07_Z_F","Rep07_Z_M",
  "Rep08_W_F","Rep08_W_M","Rep08_Z_F","Rep08_Z_M",
  "Rep09_W_F","Rep09_W_M","Rep09_Z_F","Rep09_Z_M",
  "Rep10_W_F","Rep10_W_M","Rep10_Z_F","Rep10_Z_M",
  "Rep11_W_F","Rep11_W_M","Rep11_Z_F","Rep11_Z_M",
  "Rep12_W_F","Rep12_W_M","Rep12_Z_F","Rep12_Z_M",
  "Rep13_W_F","Rep13_W_M","Rep13_Z_F","Rep13_Z_M",
  "Rep14_W_F","Rep14_W_M","Rep14_Z_F","Rep14_Z_M",
  "Rep15_W_F","Rep15_W_M","Rep15_Z_F","Rep15_Z_M"
)
