#!/usr/bin/env Rscript

# Compare adaptive vs fixed per-sample hap estimates and error diagonals
# Input: testing_positions_comparison.rds (in current working directory)
# Output:
#  - results_hap_err_diffs.csv
#  - plot_hap_diff_by_pos.png
#  - plot_err_diff_by_pos.png

suppressPackageStartupMessages({
  library(tidyverse)
})

data_file <- "testing_positions_comparison.rds"

if (!file.exists(data_file)) {
  stop("Data file not found: ", normalizePath(data_file))
}

cat("Loading data from:", data_file, "\n")
xx <- readRDS(data_file)

# Expect columns: CHROM, pos, sample(list chr), Groups(list), Haps(list), Err(list), Names(list), method
required_cols <- c("CHROM", "pos", "sample", "Haps", "Err", "method")
if (!all(required_cols %in% names(xx))) {
  stop("Input is missing required columns. Found: ", paste(names(xx), collapse = ", "))
}

# Unnest the parallel list-columns to per-sample rows
# tidyr::unnest with multiple columns keeps them aligned per row if lengths match
per_sample <- xx %>%
  dplyr::select(CHROM, pos, method, sample, Haps, Err) %>%
  tidyr::unnest(c(sample, Haps, Err)) %>%
  dplyr::rename(sample_id = sample, hap_vec = Haps, err_mat = Err)

# Sanity checks
stopifnot(is.list(per_sample$hap_vec))
stopifnot(is.list(per_sample$err_mat))

# Pivot wider to get adapt/fixed side-by-side for each (pos, sample)
pairs <- per_sample %>%
  dplyr::select(CHROM, pos, sample_id, method, hap_vec, err_mat) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = c(hap_vec, err_mat),
    names_sep = "_"
  )

# Compute differences
results <- pairs %>%
  dplyr::mutate(
    hap_diff = purrr::pmap_dbl(
      list(hap_vec_adapt, hap_vec_fixed),
      function(a, b) {
        if (is.null(a) || is.null(b)) return(NA_real_)
        a <- as.numeric(a)
        b <- as.numeric(b)
        sum(abs(a - b))
      }
    ),
    Err_diff = purrr::pmap_dbl(
      list(err_mat_adapt, err_mat_fixed),
      function(A, B) {
        if (is.null(A) || is.null(B)) return(NA_real_)
        sum(diag(as.matrix(A))) - sum(diag(as.matrix(B)))
      }
    )
  ) %>%
  dplyr::select(pos, sample = sample_id, hap_diff, Err_diff) %>%
  dplyr::arrange(pos, sample)

readr::write_csv(results, "results_hap_err_diffs.csv")
cat("Wrote results to results_hap_err_diffs.csv\n")

# Plots
colorblind_friendly_8 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#990099")

p1 <- results %>%
  ggplot2::ggplot(ggplot2::aes(x = pos, y = hap_diff, color = sample)) +
  ggplot2::geom_point(alpha = 0.9, size = 2) +
  ggplot2::labs(title = "Sum abs diff of hap estimates (adapt vs fixed)", x = "Position", y = "hap_diff") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave("plot_hap_diff_by_pos.png", p1, width = 10, height = 6, dpi = 300)
cat("Saved plot_hap_diff_by_pos.png\n")

p2 <- results %>%
  ggplot2::ggplot(ggplot2::aes(x = pos, y = Err_diff, color = sample)) +
  ggplot2::geom_point(alpha = 0.9, size = 2) +
  ggplot2::labs(title = "Signed diff of error-diagonal sums (adapt - fixed)", x = "Position", y = "Err_diff") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave("plot_err_diff_by_pos.png", p2, width = 10, height = 6, dpi = 300)
cat("Saved plot_err_diff_by_pos.png\n")

cat("Done.\n")


