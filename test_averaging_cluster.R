#!/usr/bin/env Rscript
library(tidyverse)

# Load the data
data <- readRDS('process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_chr2R.RDS')

# Get first sample data
first_sample <- data %>% 
  filter(sample == "Rep01_W_F") %>%
  arrange(pos) %>%
  slice(1:21)

cat("First 21 positions for Rep01_W_F:\n")
print(first_sample$pos)

cat("\nGroups for each position:\n")
for(i in 1:21) {
  cat("Position", first_sample$pos[i], "Groups:", first_sample$Groups[[i]], "\n")
}

cat("\nHaps for each position:\n")
for(i in 1:21) {
  cat("Position", first_sample$pos[i], "Haps:", first_sample$Haps[[i]], "\n")
}

# Check which positions have 8 distinguishable groups
has_8_groups <- map_lgl(first_sample$Groups, function(groups) {
  length(unique(groups)) == 8 && all(sort(unique(groups)) == 1:8)
})

cat("\nPositions with 8 distinguishable groups:\n")
print(which(has_8_groups))

# Get valid haps (positions with 8 groups)
valid_haps <- first_sample$Haps[has_8_groups]
cat("\nValid haps (positions with 8 groups):\n")
for(i in seq_along(valid_haps)) {
  cat("Valid position", i, "Haps:", valid_haps[[i]], "\n")
}

# Calculate average haps
if(length(valid_haps) > 0) {
  cat("\nCalculating average haps...\n")
  avg_haps <- reduce(valid_haps, `+`) / length(valid_haps)
  cat("Average haps:", avg_haps, "\n")
  cat("Sum of average haps:", sum(avg_haps), "\n")
  
  # Normalize
  avg_haps_norm <- avg_haps / sum(avg_haps)
  cat("Normalized average haps:", avg_haps_norm, "\n")
  cat("Sum of normalized haps:", sum(avg_haps_norm), "\n")
} else {
  cat("\nNo valid positions found!\n")
}
