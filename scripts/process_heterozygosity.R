#!/usr/bin/env Rscript

library(tidyverse)

#' Add genetic map positions (cM) to a dataframe
#' @param df Dataframe with chr and pos columns
#' @return Dataframe with added cM column
add_genetic <- function(df) {
  df$cM <- rep(NA, nrow(df))
  fm <- read.table("/dfs7/adl/tdlong/fly_pool/zincClean/helperfiles/flymap.r6.txt", header = FALSE)
  colnames(fm) <- c("chr", "pos", "cM")
  library(splines)
  
  for (chrs in c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")) {
    fmX <- fm %>% filter(chr == chrs)
    out <- ksmooth(fmX$pos, fmX$cM, kernel = "normal", bandwidth = 3e6)
    f_of_x <- splinefun(out$x, out$y)
    temp <- f_of_x(df$pos[df$chr == chrs])
    df$cM[df$chr == chrs] <- temp
  }
  df
}

#' Process heterozygosity files and create sliding window summaries
#' @param input_dir Directory containing nHet.*.txt files
#' @param output_file Output file name for the processed table
#' @return Processed heterozygosity table
process_heterozygosity <- function(input_dir, output_file = NULL) {
  # Ensure input directory ends with "/"
  if (!grepl("/$", input_dir)) {
    input_dir <- paste0(input_dir, "/")
  }
  
  # Set default output file if not specified
  if (is.null(output_file)) {
    output_file <- paste0(input_dir, "het.table.R")
  }
  
  ll <- list()
  for (CHR in c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")) {
    # Read heterozygosity file
    xx <- as_tibble(read.table(paste0(input_dir, "nHet.", CHR, ".txt"), header = TRUE))
    colnames(xx)[1:2] <- c("chr", "pos")
    
    # Add genetic map positions
    xx <- add_genetic(xx)
    
    # Calculate sliding window bins
    Max_cM <- max(xx$cM)
    Min_cM <- min(xx$cM)
    xx <- xx %>% 
      mutate(bincM = round(10 * (cM - Min_cM) / (Max_cM - Min_cM), 0)) %>%
      select(-cM, -chr, -pos) %>%
      pivot_longer(!bincM, names_to = "sample", values_to = "nHet") %>%
      group_by(bincM, sample) %>%
      summarize(avg_nHet = mean(nHet, na.rm = TRUE), .groups = "drop")
    
    ll[[CHR]] <- xx
  }
  
  # Combine all chromosomes
  df <- bind_rows(ll, .id = "CHR")
  
  # Write output
  write.table(df, output_file, row.names = FALSE, quote = FALSE)
  
  cat("Heterozygosity processing complete.\n")
  cat("Output written to:", output_file, "\n")
  
  return(df)
}

# Command line execution
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript process_heterozygosity.R <input_dir> [output_file]\n")
    cat("  input_dir: Directory containing nHet.*.txt files\n")
    cat("  output_file: Optional output file name (default: input_dir/het.table.R)\n")
    quit(status = 1)
  }
  
  input_dir <- args[1]
  output_file <- if (length(args) > 1) args[2] else NULL
  
  cat("Processing heterozygosity files from:", input_dir, "\n")
  het_table <- process_heterozygosity(input_dir, output_file)
  cat("Done!\n")
}
