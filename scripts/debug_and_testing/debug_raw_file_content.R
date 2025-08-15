#!/usr/bin/env Rscript

# Debug Raw File Content - Cluster Version
# This script examines the actual raw content of the REFALT file
# to understand its format and structure

cat("=== DEBUG RAW FILE CONTENT ===\n\n")

refalt_file <- "process/JUICE/RefAlt.chr2R.txt"

# 1. Check file size and basic info
cat("1. File information:\n")
file_info <- file.info(refalt_file)
cat("File size:", round(file_info$size/1e6, 2), "MB\n")
cat("File permissions:", file_info$mode, "\n")

# 2. Read first few lines as raw text
cat("\n2. First 10 lines as raw text:\n")
cat("----------------------------------------\n")
con <- file(refalt_file, "r")
for (i in 1:10) {
  line <- readLines(con, n=1)
  if (length(line) == 0) break
  cat(sprintf("Line %d (%d chars): %s\n", i, nchar(line), line))
}
close(con)

# 3. Count delimiters in first few lines
cat("\n3. Delimiter analysis:\n")
cat("----------------------------------------\n")
con <- file(refalt_file, "r")
for (i in 1:5) {
  line <- readLines(con, n=1)
  if (length(line) == 0) break
  
  # Count different delimiters
  tab_count <- lengths(regmatches(line, gregexpr("\t", line)))
  space_count <- lengths(regmatches(line, gregexpr(" ", line)))
  comma_count <- lengths(regmatches(line, gregexpr(",", line)))
  
  cat(sprintf("Line %d: %d tabs, %d spaces, %d commas\n", 
              i, tab_count, space_count, comma_count))
  
  # Show first 100 characters with visible delimiters
  line_visible <- gsub("\t", "→", line)
  line_visible <- gsub(" ", "·", line_visible)
  cat("  First 100 chars: ", substr(line_visible, 1, 100), "\n")
}
close(con)

# 4. Try different read methods
cat("\n4. Testing different read methods:\n")
cat("----------------------------------------\n")

# Method 1: readLines with different separators
cat("Method 1: readLines with tab split\n")
con <- file(refalt_file, "r")
line1 <- readLines(con, n=1)
close(con)
if (length(line1) > 0) {
  parts <- strsplit(line1, "\t")[[1]]
  cat("Line 1 split by tab:", length(parts), "parts\n")
  cat("First 5 parts:", paste(head(parts, 5), collapse=" | "), "\n")
}

# Method 2: read.table with different separators
cat("\nMethod 2: read.table with tab separator\n")
tryCatch({
  test_data <- read.table(refalt_file, nrows=1, sep="\t", fill=TRUE)
  cat("✓ Tab separator worked, columns:", ncol(test_data), "\n")
  cat("Column names:", paste(names(test_data), collapse=", "), "\n")
}, error = function(e) {
  cat("✗ Tab separator failed:", e$message, "\n")
})

# Method 3: read.table with space separator
cat("\nMethod 3: read.table with space separator\n")
tryCatch({
  test_data <- read.table(refalt_file, nrows=1, sep=" ", fill=TRUE)
  cat("✓ Space separator worked, columns:", ncol(test_data), "\n")
  cat("Column names:", paste(names(test_data), collapse=", "), "\n")
}, error = function(e) {
  cat("✗ Space separator failed:", e$message, "\n")
})

# Method 4: read.csv with comma separator
cat("\nMethod 4: read.csv with comma separator\n")
tryCatch({
  test_data <- read.csv(refalt_file, nrows=1, fill=TRUE)
  cat("✓ Comma separator worked, columns:", ncol(test_data), "\n")
  cat("Column names:", paste(names(test_data), collapse=", "), "\n")
}, error = function(e) {
  cat("✗ Comma separator failed:", e$message, "\n")
})

cat("\n=== RAW CONTENT DEBUG COMPLETE ===\n")
