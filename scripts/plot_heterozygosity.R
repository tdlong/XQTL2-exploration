#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

#' Plot heterozygosity across chromosomes with automatic founder exclusion
#' @param het_table_file Path to the heterozygosity table file
#' @param exclude_patterns Patterns to exclude (default: A and B founders)
#' @param output_file Optional output file for the plot
#' @return ggplot object
plot_heterozygosity <- function(het_table_file, 
                               exclude_patterns = c("^A[0-9]+", "^B[0-9]+", "^AB[0-9]+"),
                               output_file = NULL) {
  
  # Read the heterozygosity table
  df <- as_tibble(read.table(het_table_file, header = TRUE)) %>% 
    mutate(sample = as.factor(sample))
  
  # Filter out founder strains based on patterns
  df_filtered <- df
  for (pattern in exclude_patterns) {
    df_filtered <- df_filtered %>% 
      filter(!grepl(pattern, sample, ignore.case = TRUE))
  }
  
  # Get list of remaining samples for reference
  remaining_samples <- levels(df_filtered$sample)
  cat("Samples included in analysis:\n")
  cat(paste(remaining_samples, collapse = ", "), "\n")
  
  # Create the plot
  p <- ggplot(data = df_filtered, 
              aes(x = bincM, y = avg_nHet, group = sample, color = sample)) +
    geom_line(linewidth = 1) + 
    facet_wrap(~ CHR, ncol = 2, scales = "free_x") +
    labs(
      title = "Sliding Window Heterozygosity Analysis",
      x = "Genetic Position (cM bins)",
      y = "Average Heterozygosity",
      color = "Sample"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(6, 6, 6, 6),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "lightblue"),
      strip.text = element_text(face = "bold")
    ) +
    guides(color = guide_legend(
      override.aes = list(linewidth = 1.5),
      nrow = 2
    ))
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 12, height = 8, dpi = 300)
    cat("Plot saved to:", output_file, "\n")
  }
  
  return(p)
}

#' Quick heterozygosity analysis with common strain groupings
#' @param het_table_file Path to the heterozygosity table file
#' @param strain_groups Named list of strain groups to analyze
#' @param output_dir Directory to save plots
analyze_strain_groups <- function(het_table_file, 
                                 strain_groups = NULL,
                                 output_dir = ".") {
  
  # Default strain groups if none provided
  if (is.null(strain_groups)) {
    # You can customize these default groups based on your common analyses
    strain_groups <- list(
      "test_strains" = c("strain1", "strain2", "strain3"),
      "control_strains" = c("control1", "control2")
    )
  }
  
  # Create plots for each group
  for (group_name in names(strain_groups)) {
    strains <- strain_groups[[group_name]]
    
    # Filter data for this group
    df <- as_tibble(read.table(het_table_file, header = TRUE)) %>% 
      mutate(sample = as.factor(sample)) %>%
      filter(sample %in% strains)
    
    if (nrow(df) > 0) {
      # Create plot for this group
      p <- ggplot(data = df, 
                  aes(x = bincM, y = avg_nHet, group = sample, color = sample)) +
        geom_line(linewidth = 1) + 
        facet_wrap(~ CHR, ncol = 2, scales = "free_x") +
        labs(
          title = paste("Heterozygosity Analysis:", group_name),
          x = "Genetic Position (cM bins)",
          y = "Average Heterozygosity",
          color = "Sample"
        ) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          panel.grid.minor = element_blank()
        )
      
      # Save plot
      output_file <- file.path(output_dir, paste0("heterozygosity_", group_name, ".png"))
      ggsave(output_file, p, width = 12, height = 8, dpi = 300)
      cat("Plot saved for", group_name, "to:", output_file, "\n")
    }
  }
}

# Command line execution
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript plot_heterozygosity.R <het_table_file> [output_file]\n")
    cat("  het_table_file: Path to the heterozygosity table (e.g., het.table.R)\n")
    cat("  output_file: Optional output file for the plot (default: heterozygosity_plot.png)\n")
    cat("\nExample:\n")
    cat("  Rscript plot_heterozygosity.R process/heterozygosity_analysis/het.table.R\n")
    cat("  Rscript plot_heterozygosity.R het.table.R my_plot.png\n")
    quit(status = 1)
  }
  
  het_table_file <- args[1]
  output_file <- if (length(args) > 1) args[2] else "heterozygosity_plot.png"
  
  cat("Creating heterozygosity plot from:", het_table_file, "\n")
  p <- plot_heterozygosity(het_table_file, output_file = output_file)
  cat("Plot created successfully!\n")
}
