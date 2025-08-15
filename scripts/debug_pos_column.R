#!/usr/bin/env Rscript

# Debug POS Column - Simple inspection of the POS column data
cat("=== DEBUG POS COLUMN ===\n\n")

filein <- "process/JUICE/RefAlt.chr2R.txt"

# Read first few lines to see raw content
cat("1. First 5 lines of raw file:\n")
con <- file(filein, "r")
for (i in 1:5) {
  line <- readLines(con, n=1)
  cat(sprintf("Line %d: %s\n", i, line))
}
close(con)

# Try to read with different approaches
cat("\n2. Testing different read approaches:\n")

# Approach 1: read.table with tab separator
cat("\nApproach 1: read.table with tab separator\n")
tryCatch({
  df1 <- read.table(filein, nrows=5, sep="\t", header=TRUE)
  cat("✓ Success - columns:", ncol(df1), "\n")
  cat("Column names:", paste(names(df1), collapse=", "), "\n")
  cat("POS column values:", paste(df1$POS, collapse=", "), "\n")
  cat("POS column class:", class(df1$POS), "\n")
}, error = function(e) {
  cat("✗ Failed:", e$message, "\n")
})

# Approach 2: read.table with space separator
cat("\nApproach 2: read.table with space separator\n")
tryCatch({
  df2 <- read.table(filein, nrows=5, sep=" ", header=TRUE)
  cat("✓ Success - columns:", ncol(df2), "\n")
  cat("Column names:", paste(names(df2), collapse=", "), "\n")
  if ("POS" %in% names(df2)) {
    cat("POS column values:", paste(df2$POS, collapse=", "), "\n")
    cat("POS column class:", class(df2$POS), "\n")
  } else {
    cat("No POS column found\n")
  }
}, error = function(e) {
  cat("✗ Failed:", e$message, "\n")
})

# Approach 3: read.csv with comma separator
cat("\nApproach 3: read.csv with comma separator\n")
tryCatch({
  df3 <- read.csv(filein, nrows=5, header=TRUE)
  cat("✓ Success - columns:", ncol(df3), "\n")
  cat("Column names:", paste(names(df3), collapse=", "), "\n")
  if ("POS" %in% names(df3)) {
    cat("POS column values:", paste(df3$POS, collapse=", "), "\n")
    cat("POS column class:", class(df3$POS), "\n")
  } else {
    cat("No POS column found\n")
  }
}, error = function(e) {
  cat("✗ Failed:", e$message, "\n")
})

cat("\n=== POS COLUMN DEBUG COMPLETE ===\n")
