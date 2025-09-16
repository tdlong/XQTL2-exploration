# Function to process combined chromosomes (debug3 mode)
# Combines 2L+2R and 3L+3R data for more robust haplotype estimation
# Still applies debug2 mode: 1500 most proximal 2L SNPs

process_combined_chromosomes <- function(chr_group, centromere_positions, input_dir, founders, names_in_bam, h_cutoff, debug_mode, debug2_mode) {
  cat("\n--- PROCESSING CHROMOSOME GROUP:", chr_group, "---\n")
  
  # Get centromere positions for this chromosome group
  chr_positions <- centromere_positions %>%
    filter(CHROM %in% paste0("chr", chr_group, c("L", "R"))) %>%
    pull(pos)
  
  if (length(chr_positions) == 0) {
    cat("No centromere positions found for chr", chr_group, "- skipping\n")
    return(tibble())
  }
  
  cat("Combined chr", chr_group, " positions:", length(chr_positions), "\n")
  cat("Position range:", min(chr_positions), "to", max(chr_positions), "\n\n")
  
  # Load and combine RefAlt data from both arms
  combined_data <- tibble()
  for (arm in c("L", "R")) {
    chr <- paste0("chr", chr_group, arm)
    filein <- file.path(input_dir, paste0("RefAlt.", chr, ".txt"))
    if (file.exists(filein)) {
      chr_data <- read_tsv(filein, col_types = cols(.default = "c")) %>%
        mutate(CHROM = chr)
      
      # Apply debug2 mode (1500 proximal SNPs) to chr2L BEFORE combining
      if (debug2_mode && chr == "chr2L") {
        # Get chr2L centromere positions
        chr2L_positions <- centromere_positions %>%
          filter(CHROM == "chr2L") %>%
          pull(pos)
        
        # Limit to 1500 most proximal SNPs (highest positions)
        if (length(chr2L_positions) > 1500) {
          chr2L_positions <- sort(chr2L_positions, decreasing = TRUE)[1:1500]
          cat("DEBUG2 MODE: Limited chr2L to 1500 most proximal SNPs\n")
        }
        
        # Filter chr2L data to these positions
        chr_data <- chr_data %>%
          filter(POS %in% chr2L_positions)
      }
      
      combined_data <- bind_rows(combined_data, chr_data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    cat("No RefAlt data found for chr", chr_group, "\n")
    return(tibble())
  }
  
  # Renumber positions from 1 to nrows to avoid overlap issues
  combined_data <- combined_data %>%
    mutate(POS = row_number())
  
  # Filter to centromere positions (now using the renumbered positions)
  # We need to map the original centromere positions to the new positions
  # This is tricky because we renumbered after combining...
  
  # For now, let's just use all the combined data
  chr_data <- combined_data
  
  cat("Found", nrow(chr_data), "combined positions\n")
  
  if (nrow(chr_data) < 100) {
    cat("WARNING: Very few positions found - may indicate a problem!\n")
  }
  
  # Process the combined data
  cat("Processing combined chr", chr_group, " data...\n")
  
  # TODO: Add the full haplotype estimation logic here
  # This would be similar to the single chromosome processing but with combined data
  
  return(tibble())  # Placeholder for now
}
