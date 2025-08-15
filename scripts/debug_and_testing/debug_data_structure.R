#!/usr/bin/env Rscript

# Debug Data Structure - Cluster Version
# This script checks the actual structure of REFALT data on the cluster
# to identify why founder columns are missing

cat("=== DEBUG DATA STRUCTURE ON CLUSTER ===\n\n")

# 1. Check if parameter file exists and load it
cat("1. Loading parameter file...\n")
param_file <- "helpfiles/JUICE/JUICE_haplotype_parameters.R"
if (file.exists(param_file)) {
  cat("✓ Parameter file found:", param_file, "\n")
  source(param_file)
  cat("✓ Founders loaded:", paste(founders, collapse=", "), "\n")
} else {
  cat("✗ Parameter file NOT found:", param_file, "\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Files in helpfiles/:\n")
  list.files("helpfiles/", recursive=TRUE)
  quit(status=1)
}

# 2. Check if REFALT file exists
cat("\n2. Checking REFALT file...\n")
refalt_file <- "process/JUICE/RefAlt.chr2R.txt"
if (file.exists(refalt_file)) {
  cat("✓ REFALT file found:", refalt_file, "\n")
  file_info <- file.info(refalt_file)
  cat("File size:", round(file_info$size/1e6, 2), "MB\n")
} else {
  cat("✗ REFALT file NOT found:", refalt_file, "\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Files in process/:\n")
  list.files("process/", recursive=TRUE)
  quit(status=1)
}

# 3. Load a small sample of data to check structure
cat("\n3. Loading data sample...\n")
tryCatch({
  # Load first 100 rows to see structure
  data_sample <- read.table(refalt_file, nrows=100, header=TRUE, sep="\t")
  cat("✓ Data loaded successfully\n")
  cat("Data dimensions:", dim(data_sample), "\n")
  cat("Data columns:", paste(names(data_sample), collapse=", "), "\n")
  
  # Check if founder columns exist
  cat("\n4. Checking founder columns...\n")
  missing_founders <- founders[!founders %in% names(data_sample)]
  if (length(missing_founders) == 0) {
    cat("✓ All founder columns found!\n")
  } else {
    cat("✗ Missing founder columns:", paste(missing_founders, collapse=", "), "\n")
  }
  
  # Show first few rows
  cat("\n5. First few rows of data:\n")
  print(head(data_sample[, 1:min(10, ncol(data_sample))]))
  
  # Check data types
  cat("\n6. Column data types:\n")
  col_types <- sapply(data_sample, class)
  print(col_types)
  
}, error = function(e) {
  cat("✗ Error loading data:", e$message, "\n")
  cat("Trying to read without header...\n")
  
  # Try without header
  tryCatch({
    data_sample <- read.table(refalt_file, nrows=100, header=FALSE, sep="\t")
    cat("✓ Data loaded without header\n")
    cat("Data dimensions:", dim(data_sample), "\n")
    cat("First few columns (showing first 10):\n")
    print(head(data_sample[, 1:min(10, ncol(data_sample))]))
  }, error = function(e2) {
    cat("✗ Still failed:", e2$message, "\n")
  })
})

cat("\n=== DEBUG COMPLETE ===\n")
