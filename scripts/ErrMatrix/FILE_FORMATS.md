# File Format Documentation

## Original Adaptive Results (`adaptive_window_h4_results_chr3R.RDS`)

**Format**: Nested tibble with one row per (CHROM, pos, sample) combination, but data stored in list columns

**Structure**:
- **Rows**: 163,740 (one per sample per position)
- **Columns**: 7
  - `CHROM`: character - chromosome name (e.g., "chr3R")
  - `pos`: numeric - genomic position
  - `sample`: character - sample identifier (e.g., "Rep01_W_F")
  - `Groups`: list with 1 element - group assignments (numeric vector)
  - `Haps`: list with 1 element - haplotype frequency estimates (numeric vector, length = 8)
  - `Err`: list with 1 element - 8x8 error/covariance matrix for haplotype estimates
  - `Names`: list with 1 element - founder names (character vector, length = 8)

**Key Properties**:
- Data is stored in list columns (need unnesting to access)
- `Err[[1]]` matrices have proper row/column names matching `Names[[1]]`
- `Haps[[1]]` vectors have length 8 (one per founder)
- `Names[[1]]` vectors have length 8 with founder identifiers

## Reshaped Adaptive Results (`R.haps.chr3R.out.rds`)

**Format**: Nested tibble with one row per (CHROM, pos) combination, samples stored as lists

**Structure**:
- **Rows**: 2,729 (one per position)
- **Columns**: 7
  - `CHROM`: character - chromosome name (e.g., "chr3R")
  - `pos`: numeric - genomic position
  - `sample`: list of character vectors - sample identifiers (length = 60, one per sample)
  - `Groups`: list of numeric vectors - group assignments (length = 60, one per sample)
  - `Haps`: list of numeric vectors - haplotype frequency estimates (length = 60, each vector has length 8)
  - `Err`: list of numeric matrices - error/covariance matrices (length = 60, each matrix is 8x8)
  - `Names`: list of character vectors - founder names (length = 60, each vector has length 8)

**Key Properties**:
- Nested format - each row contains lists of 60 samples
- `Err` matrices have proper row/column names matching `Names`
- `Haps` vectors have length 8 (one per founder)
- `Names` vectors have length 8 with founder identifiers
- **Requires unnesting** to get one row per sample for comparison

## Comparison Requirements

To compare these files:

1. **Original file**: Must be unnested using `tidyr::unnest(c(sample, Groups, Haps, Err, Names))`
2. **Reshaped file**: Must be unnested using `tidyr::unnest(c(sample, Groups, Haps, Err, Names))`
3. **After unnesting**: Both files have identical structure with one row per (CHROM, pos, sample)
4. **Data integrity**: The unnested reshaped data should be identical to the unnested original data

## Data Flow

```
Original Adaptive Results (163,740 rows)
    ↓
Reshaping (group_by CHROM, pos, summarise with lists)
    ↓
Reshaped Results (2,729 rows with list columns)
    ↓
Unnesting (tidyr::unnest)
    ↓
Unnested Reshaped Results (163,740 rows)
    ↓
Should be identical to Original Adaptive Results
```

## File Creation Code

**Original file creation** (in `BASE_VAR_WIDE.R`):
```r
# Results are already in unnested format from est_haps_var()
adaptive_results <- est_haps_var(df4, founders, ...)
saveRDS(adaptive_results, adaptive_original_file)
```

**Reshaped file creation** (in `BASE_VAR_WIDE.R`):
```r
# Reshape to one row per position with samples as lists
adaptive_data_reshaped <- adaptive_results %>%
  dplyr::arrange(CHROM, pos, sample) %>%
  dplyr::group_by(CHROM, pos) %>%
  dplyr::summarise({
    ord <- order(sample)
    tibble(
      sample = list(sample[ord]),
      Groups = list(Groups[ord]),
      Haps   = list(Haps[ord]),
      Err    = list(Err[ord]),
      Names  = list(Names[ord])
    )
  }, .groups = "drop")
saveRDS(adaptive_data_reshaped, adaptive_reshaped_file)
```

## Validation

The reshaped file should preserve:
- All data values (no recomputation)
- Sample ordering within each position
- Matrix dimensions and names
- Vector lengths and content
- Founder name alignment between `Names` and `Err` matrix dimensions