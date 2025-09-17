#!/usr/bin/env Rscript

# Minimal extractor: save ONLY the inputs needed by the haplotype estimator
# No estimation, no clustering. You can download the RDS and run locally.
#
# Usage:
# Rscript scripts/ErrMatrix/extract_estimator_inputs.R <chr> <param_file> <output_dir> <position> <sample_name> <window_size_bp>
# Example:
# Rscript scripts/ErrMatrix/extract_estimator_inputs.R chr3R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 19780000 Rep01_W_F 100000

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript scripts/ErrMatrix/extract_estimator_inputs.R <chr> <param_file> <output_dir> <position> <sample_name> <window_size_bp>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
test_position <- as.numeric(args[4])
sample_name <- args[5]
window_size <- as.numeric(args[6])

cat("=== EXTRACT ESTIMATOR INPUTS ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output dir:", output_dir, "\n")
cat("Position:", test_position, "\n")
cat("Sample:", sample_name, "\n")
cat("Window size:", window_size, "bp\n\n")

# Load parameters (expects 'founders' and 'h_cutoff' at minimum)
source(param_file)
if (!exists("founders")) stop("'founders' not found in param file")

# Read REF/ALT counts for chromosome
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) stop("RefAlt file not found: ", refalt_file)

cat("Loading REFALT from:", refalt_file, "\n")
refalt_data <- readr::read_tsv(
  refalt_file,
  col_names = c("chr", "pos", "name", "ref_count", "alt_count"),
  col_types = "cdcdd",
  show_col_types = FALSE
)

# Convert to frequencies (ALT / total) and keep founders + target sample only within window
window_start <- test_position - window_size / 2
window_end <- test_position + window_size / 2

observed_data <- refalt_data %>%
  dplyr::mutate(total_count = .data$ref_count + .data$alt_count,
                freq = .data$alt_count / .data$total_count) %>%
  dplyr::filter(.data$total_count > 0) %>%
  dplyr::select(.data$pos, .data$name, .data$freq) %>%
  dplyr::filter(.data$pos >= window_start & .data$pos <= window_end &
                .data$name %in% c(founders, sample_name))

cat("SNPs in window:", nrow(observed_data), "\n")
if (nrow(observed_data) == 0) stop("No SNPs in window.")

# Wide matrix: rows=SNPs, columns=founders + sample
wide_data <- observed_data %>%
  tidyr::pivot_wider(names_from = .data$name, values_from = .data$freq) %>%
  dplyr::arrange(.data$pos)

missing_cols <- setdiff(c(founders, sample_name), names(wide_data))
if (length(missing_cols) > 0) {
  stop("Missing expected columns in wide data: ", paste(missing_cols, collapse = ", "))
}

founder_matrix <- wide_data %>% dplyr::select(dplyr::all_of(founders)) %>% as.matrix()
sample_freqs <- wide_data %>% dplyr::pull(!!sample_name)

# Clean rows with any NA in founders or sample
complete_rows <- stats::complete.cases(founder_matrix) & !is.na(sample_freqs)
founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
sample_freqs_clean <- sample_freqs[complete_rows]
positions_clean <- wide_data$pos[complete_rows]

cat("Clean matrix:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
if (nrow(founder_matrix_clean) < 1) stop("No complete rows after cleaning.")

# Save minimal inputs
out <- list(
  chr = chr,
  position = test_position,
  window_size_bp = window_size,
  founders = founders,
  sample_name = sample_name,
  positions = positions_clean,
  founder_matrix = founder_matrix_clean,
  sample_freqs = sample_freqs_clean
)

outfile <- sprintf("estimator_inputs_%s_%d_%s.win%d.rds", chr, test_position, sample_name, window_size)
saveRDS(out, outfile)
cat("Saved:", outfile, "\n")


