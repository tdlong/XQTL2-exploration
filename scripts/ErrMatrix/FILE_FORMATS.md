# File Format Documentation

## Adaptive Window Results (Original Format)
**File**: `process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_chr3R.RDS`

**Structure**: One row per (CHROM, pos, sample)
- `CHROM`: character (e.g., "chr3R")
- `pos`: numeric (genomic position)
- `sample`: character (e.g., "Rep01_W_F")
- `Groups`: list of integers (e.g., `c(1,2,3,4,5,6,7,8)`)
- `Haps`: list of named numeric vectors (e.g., `c(B1=0.1, B2=0.2, ...)`)
- `Err`: list of 8x8 matrices with row/column names matching founder names
- `Names`: list of character vectors (e.g., `c("B1","B2","B3","B4","B5","B6","B7","AB8")`)

**Key Properties**:
- `Err[[i]]` is a matrix with `rownames(Err[[i]]) == Names[[i]]` and `colnames(Err[[i]]) == Names[[i]]`
- `Haps[[i]]` is a named vector with `names(Haps[[i]]) == Names[[i]]`
- `Groups[[i]]` indicates which founders are grouped together (1-8 for full separation)

## Reshaped Adapt H4 Results
**File**: `process/ZINC2/haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds`

**Structure**: One row per (CHROM, pos) with nested lists for samples
- `CHROM`: character (e.g., "chr3R")
- `pos`: numeric (genomic position)
- `sample`: list of character vectors (e.g., `list(c("Rep01_W_F", "Rep01_W_M", ...))`)
- `Groups`: list of lists (e.g., `list(list(c(1,2,3,4,5,6,7,8), c(1,2,3,4,5,6,7,8), ...))`)
- `Haps`: list of lists (e.g., `list(list(c(B1=0.1, B2=0.2, ...), c(B1=0.1, B2=0.2, ...), ...))`)
- `Err`: list of lists (e.g., `list(list(matrix(...), matrix(...), ...))`)
- `Names`: list of lists (e.g., `list(list(c("B1","B2",...), c("B1","B2",...), ...))`)

**After unnesting** (via `tidyr::unnest(c(sample, Groups, Haps, Err, Names))`):
- `sample`: character (e.g., "Rep01_W_F")
- `Groups`: list of integers (e.g., `c(1,2,3,4,5,6,7,8)`)
- `Haps`: list of named numeric vectors (e.g., `c(B1=0.1, B2=0.2, ...)`)
- `Err`: list of 8x8 matrices
- `Names`: list of character vectors (e.g., `c("B1","B2","B3","B4","B5","B6","B7","AB8")`)

**Key Properties**:
- After unnesting, structure should match original format
- `Err[[i]]` should be a matrix with proper row/column names
- `Haps[[i]]` should be a named vector
- `Names[[i]]` should be a character vector

## Expected Alignment
For comparison, we expect:
1. `Names_o[[i]] == Names_r[[i]]` (same founder order)
2. `Err_o[[i]]` and `Err_r[[i]]` are both 8x8 matrices
3. `rownames(Err_r[[i]]) == Names_r[[i]]` and `colnames(Err_r[[i]]) == Names_r[[i]]`
4. `names(Haps_r[[i]]) == Names_r[[i]]`

## Current Issue
The comparison script shows:
- `has_names: TRUE` - Names are present and match
- `has_err_mats: FALSE` - Error matrices are missing or malformed
- `names_in_Er: FALSE` - The reshaped error matrices don't have expected row/column names

This suggests the reshaped `Err` matrices are either:
1. Not matrices at all (maybe lists or vectors)
2. Missing row/column names entirely
3. Have different row/column names than expected
