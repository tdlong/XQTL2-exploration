#!/usr/bin/env Rscript

# Compare old (standard) format vs new (list) format haplotype results
# Validates that haplotype frequencies (B1..AB8) are consistent

suppressPackageStartupMessages({
	library(tidyverse)
})

# Args: <chromosome> [position] [epsilon]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
	stop("Usage: Rscript scripts/debug/compare_old_vs_list_format.R <chromosome> [position] [epsilon]")
}
chromosome <- args[1]
position <- if (length(args) > 1 && nzchar(args[2])) suppressWarnings(as.numeric(args[2])) else NULL
epsilon <- if (length(args) > 2 && nzchar(args[3])) suppressWarnings(as.numeric(args[3])) else 1e-6
if (is.na(epsilon)) epsilon <- 1e-6

message("Comparing old vs list formats for ", chromosome, if (!is.null(position)) paste0(" at ", position) else "",
        ", eps=", format(epsilon, scientific = TRUE))

# File paths (relative)
old_adaptive_file <- paste0("process/ZINC2/haplotype_results/adaptive_window_h4_results_", chromosome, ".RDS")
old_smooth_file   <- paste0("process/ZINC2/haplotype_results/smooth_h4_results_", chromosome, ".RDS")
new_adaptive_file <- paste0("process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_", chromosome, ".RDS")
new_smooth_file   <- paste0("process/ZINC2/haplotype_results_list_format/smooth_h4_results_", chromosome, ".RDS")

# Helper: check file exists
check_file <- function(path) {
	if (!file.exists(path)) stop("Missing file: ", path)
}

# Load files
check_file(old_adaptive_file)
check_file(new_adaptive_file)
old_adaptive <- readRDS(old_adaptive_file)
new_adaptive <- readRDS(new_adaptive_file)

# Old smooth is optional; new smooth expected
has_old_smooth <- file.exists(old_smooth_file)
has_new_smooth <- file.exists(new_smooth_file)
if (has_old_smooth) old_smooth <- readRDS(old_smooth_file)
if (has_new_smooth) new_smooth <- readRDS(new_smooth_file)

# Extractors
founder_cols <- c("B1","B2","B3","B4","B5","B6","B7","AB8")

extract_old <- function(df, pos = NULL) {
	df2 <- df %>%
		{
			if (!is.null(pos)) dplyr::filter(., .data$pos == !!pos) else .
		} %>%
		dplyr::select(chr, pos, sample, dplyr::all_of(founder_cols)) %>%
		dplyr::rename(CHROM = chr) %>%
		dplyr::arrange(CHROM, pos, sample)
	return(df2)
}

extract_new <- function(df, pos = NULL) {
	df2 <- df %>%
		{
			if (!is.null(pos)) dplyr::filter(., .data$pos == !!pos) else .
		} %>%
		dplyr::select(CHROM, pos, sample, Haps) %>%
		dplyr::mutate(
			B1 = purrr::map_dbl(Haps, ~ .x[1]),
			B2 = purrr::map_dbl(Haps, ~ .x[2]),
			B3 = purrr::map_dbl(Haps, ~ .x[3]),
			B4 = purrr::map_dbl(Haps, ~ .x[4]),
			B5 = purrr::map_dbl(Haps, ~ .x[5]),
			B6 = purrr::map_dbl(Haps, ~ .x[6]),
			B7 = purrr::map_dbl(Haps, ~ .x[7]),
			AB8 = purrr::map_dbl(Haps, ~ .x[8])
		) %>%
		dplyr::select(-Haps) %>%
		dplyr::arrange(CHROM, pos, sample)
	return(df2)
}

# Comparator
compare_frames <- function(old_df, new_df, label, eps) {
	message("\n===== ", label, " =====")
	# Align rows that exist in both
	joined <- dplyr::inner_join(
		old_df, new_df,
		by = c("CHROM","pos","sample"), suffix = c("_old","_new")
	)
	message("Rows in common: ", nrow(joined))
	if (nrow(joined) == 0) {
		message("No overlapping rows; skipping.")
		return(invisible(NULL))
	}
	# Per-founder diffs
	diff_df <- purrr::map_dfc(founder_cols, function(fcol) {
		old_col <- joined[[paste0(fcol, "_old")]]
		new_col <- joined[[paste0(fcol, "_new")]]
		setNames(list(old_col - new_col), fcol)
	})
	# Only consider rows with complete data across all 8 founders
	valid_rows <- stats::complete.cases(diff_df)
	ssq <- rep(NA_real_, nrow(joined))
	ssq[valid_rows] <- rowSums(as.matrix(diff_df[valid_rows, , drop = FALSE])^2)
	n_examined <- sum(!is.na(ssq))
	n_exceed <- sum(ssq > eps, na.rm = TRUE)
	pct <- if (n_examined > 0) 100 * n_exceed / n_examined else NA_real_
	max_ssq <- if (n_examined > 0) max(ssq, na.rm = TRUE) else NA_real_
	message(sprintf("SSQ>eps: examined=%d, exceed=%d (%.4f%%), max_ssq=%.3e", n_examined, n_exceed, pct, max_ssq))
	invisible(list(n_examined = n_examined, n_exceed = n_exceed, max_ssq = max_ssq))
}

# Run comparisons
old_adaptive_freq <- extract_old(old_adaptive, position)
new_adaptive_freq <- extract_new(new_adaptive, position)
compare_frames(old_adaptive_freq, new_adaptive_freq, "ADAPTIVE_H4", epsilon)

if (has_old_smooth && has_new_smooth) {
	old_smooth_freq <- extract_old(old_smooth, position)
	new_smooth_freq <- extract_new(new_smooth, position)
	compare_frames(old_smooth_freq, new_smooth_freq, "SMOOTH_H4", epsilon)
} else {
	message("\nSMOOTH_H4: one or both files missing; skipping")
}

message("\nDone.")
