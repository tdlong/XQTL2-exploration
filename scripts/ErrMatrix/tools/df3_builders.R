#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# process_refalt_data(): copied from working production code to avoid sourcing
# the full BASE_VAR_WIDE.R. This constructs df3 in WIDE format exactly as used
# by the estimator pipeline.

process_refalt_data <- function(refalt_file, founders) {
  cat("Loading RefAlt data from:", refalt_file, "\n")
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)

  df3 <- refalt_data %>%
    pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
    mutate(
      RefAlt = stringr::str_sub(lab, 1, 3),
      name = stringr::str_sub(lab, 5)
    ) %>%
    select(-lab) %>%
    pivot_wider(names_from = RefAlt, values_from = count) %>%
    mutate(
      freq = ALT / (REF + ALT),
      total_count = REF + ALT
    ) %>%
    select(-REF, -ALT)

  founder_wide <- df3 %>%
    filter(name %in% founders) %>%
    dplyr::select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)

  quality_filtered_positions <- founder_wide %>%
    mutate(across(all_of(founders), ~ .x)) %>%
    rowwise() %>%
    mutate(
      zeros = sum(is.na(c_across(all_of(founders))) | is.infinite(c_across(all_of(founders))))
    ) %>%
    ungroup() %>%
    filter(zeros == 0) %>%
    pull(POS)

  df3 <- df3 %>%
    filter(POS %in% quality_filtered_positions)

  df3 <- df3 %>% dplyr::select(POS, name, freq) %>% pivot_wider(names_from = name, values_from = freq)

  cat("✓ Processed", nrow(df3), "rows for", ncol(df3)-1, "samples in WIDE format\n")
  return(df3)
}


