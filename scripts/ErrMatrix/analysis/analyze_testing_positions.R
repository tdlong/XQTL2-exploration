#!/usr/bin/env Rscript

# Analyze Testing Positions for Bug Detection
# This script analyzes the extracted testing positions to identify differences
# between adaptive and fixed window methods

# Load required libraries
library(tidyverse)

# Colorblind-friendly palette
colorblind_friendly_8 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#990099")

# Load the extracted data
data_file <- "testing_positions_comparison.rds"

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file, "\nPlease run the extraction script first.")
}

cat("Loading extracted data...\n")
data <- readRDS(data_file)

cat("Data loaded successfully!\n")
cat("Shape:", nrow(data), "rows,", ncol(data), "columns\n")
cat("Methods:", paste(unique(data$method), collapse = ", "), "\n")
cat("Positions:", paste(unique(data$pos), collapse = ", "), "\n\n")

# Display data structure
cat("Data structure:\n")
str(data)

# Show first few rows
cat("\nFirst 10 rows:\n")
print(head(data, 10))

# Basic summary by method and position
cat("\nSummary by method and position:\n")
summary_data <- data %>%
  dplyr::group_by(method, pos) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(pos, method)
print(summary_data)

# Check for missing positions
all_positions <- c(19780000, 19790000, 19800000, 19810000, 19820000, 19830000, 19840000)
missing_positions <- setdiff(all_positions, unique(data$pos))
if (length(missing_positions) > 0) {
  cat("\nWARNING: Missing positions:", paste(missing_positions, collapse = ", "), "\n")
}

# Check for differences in column structure between methods
adapt_cols <- names(data[data$method == "adapt", ])
fixed_cols <- names(data[data$method == "fixed", ])

if (!setequal(adapt_cols, fixed_cols)) {
  cat("\nWARNING: Different columns between methods\n")
  cat("Adaptive columns:", paste(adapt_cols, collapse = ", "), "\n")
  cat("Fixed columns:", paste(fixed_cols, collapse = ", "), "\n")
}

# Identify numeric columns for comparison
numeric_cols <- data %>%
  dplyr::select_if(is.numeric) %>%
  names() %>%
  setdiff("pos")  # Exclude position column

cat("\nNumeric columns for comparison:", paste(numeric_cols, collapse = ", "), "\n")

# Compare values between methods for each position
if (length(numeric_cols) > 0) {
  cat("\nComparing numeric values between methods:\n")
  
  for (pos in unique(data$pos)) {
    cat("\n--- Position", pos, "---\n")
    
    pos_data <- data %>%
      dplyr::filter(pos == !!pos) %>%
      dplyr::arrange(method)
    
    if (nrow(pos_data) == 2) {
      adapt_row <- pos_data[pos_data$method == "adapt", ]
      fixed_row <- pos_data[pos_data$method == "fixed", ]
      
      for (col in numeric_cols) {
        adapt_val <- adapt_row[[col]]
        fixed_val <- fixed_row[[col]]
        
        if (length(adapt_val) == 1 && length(fixed_val) == 1) {
          diff <- adapt_val - fixed_val
          cat(sprintf("%s: adapt=%.6f, fixed=%.6f, diff=%.6f\n", 
                      col, adapt_val, fixed_val, diff))
        }
      }
    } else {
      cat("Expected 2 rows (one per method), found", nrow(pos_data), "\n")
    }
  }
}

# Create comparison plots if we have enough data
if (length(numeric_cols) > 0 && length(unique(data$pos)) > 1) {
  cat("\nCreating comparison plots...\n")
  
  # Plot 1: Compare values across positions for each numeric column
  for (col in numeric_cols[1:min(4, length(numeric_cols))]) {  # Limit to first 4 columns
    p <- data %>%
      ggplot2::ggplot(aes(x = factor(pos), y = .data[[col]], color = method, shape = method)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::geom_line(aes(group = method), alpha = 0.6) +
      ggplot2::scale_color_manual(values = colorblind_friendly_8[1:2]) +
      ggplot2::scale_shape_manual(values = c(16, 17)) +
      ggplot2::labs(
        title = paste("Comparison of", col, "across testing positions"),
        x = "Position",
        y = col,
        color = "Method",
        shape = "Method"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(paste0("comparison_", col, ".png"), p, width = 10, height = 6, dpi = 300)
    cat("Saved plot: comparison_", col, ".png\n")
  }
  
  # Plot 2: Scatter plot comparing methods
  if (length(numeric_cols) >= 2) {
    col1 <- numeric_cols[1]
    col2 <- numeric_cols[2]
    
    p <- data %>%
      ggplot2::ggplot(aes(x = .data[[col1]], y = .data[[col2]], color = method, shape = method)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::geom_text(aes(label = pos), hjust = -0.2, vjust = 0.2, size = 3) +
      ggplot2::scale_color_manual(values = colorblind_friendly_8[1:2]) +
      ggplot2::scale_shape_manual(values = c(16, 17)) +
      ggplot2::labs(
        title = paste("Scatter plot:", col1, "vs", col2),
        x = col1,
        y = col2,
        color = "Method",
        shape = "Method"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    
    ggsave("scatter_comparison.png", p, width = 8, height = 6, dpi = 300)
    cat("Saved plot: scatter_comparison.png\n")
  }
}

# Save detailed comparison table
comparison_table <- data %>%
  dplyr::select(pos, method, everything()) %>%
  dplyr::arrange(pos, method)

write.csv(comparison_table, "testing_positions_detailed.csv", row.names = FALSE)
cat("\nSaved detailed comparison table: testing_positions_detailed.csv\n")

cat("\nAnalysis complete!\n")
cat("Files created:\n")
cat("- testing_positions_detailed.csv: Detailed comparison table\n")
if (length(numeric_cols) > 0) {
  cat("- comparison_*.png: Comparison plots\n")
  cat("- scatter_comparison.png: Scatter plot comparison\n")
}
