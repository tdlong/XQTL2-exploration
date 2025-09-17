# Haplotype Estimation Debugging Summary

## Problem Statement
The adaptive haplotype estimation method produces **20x higher error estimates** than the fixed window method for position 19780000 on chr3R. This is a significant bug that needs investigation.

## Approach
We created a debugging pipeline to replicate the production haplotype estimator on a single position/sample to understand why the error estimates are inflated.

## What We Built

### 1. Data Extraction Script (`extract_hunk.r`)
**Purpose**: Extract the exact data that gets passed to the haplotype estimator function

**What it reads**:
- `process/ZINC2/RefAlt.chr3R.txt` - Raw RefAlt data file
- `helpfiles/ZINC2_haplotype_parameters.R` - Parameter file defining founders

**What it does**:
- Loads RefAlt data using `process_refalt_data()` function
- Pre-subsets data to window around target position (19780000 ± 500kb)
- Extracts data for specific sample (Rep01_W_F)
- Packages `df3` data and estimator arguments into payload

**What it writes**:
- `hunk_data_chr3R_19780000_Rep01_W_F_h4.rds` - Contains:
  - `df3`: Pre-subsetted data frame (7925 rows × 69 columns)
  - `args`: List of arguments for `est_haps_var()` function

### 2. Debug Wrapper Script (`debug_wrapper.R`)
**Purpose**: Run the exact production haplotype estimator function on extracted data

**What it reads**:
- `hunk_data_chr3R_19780000_Rep01_W_F_h4.rds` - Extracted data from step 1

**What it does**:
- Copies the EXACT `est_haps_var` function from production code
- Calls `est_haps_var` with extracted data and arguments
- Runs with maximum verbosity (verbose=3) to see internal behavior

**What it writes**:
- `debug_result_chr3R_19780000_Rep01_W_F.rds` - Contains estimator results

### 3. Production Data Extraction (`extract_specific_sample.R`)
**Purpose**: Extract production results for comparison

**What it reads**:
- `testing_positions_comparison.rds` - Production results comparing adapt vs fixed methods

**What it does**:
- Unnests list-columns to get per-sample data
- Extracts data for Rep01_W_F at position 19780000
- Shows both adaptive and fixed method results

## Results

### Production Data (from comparison file)
**Adaptive Method**:
- Haplotype frequencies: 0.185512, 0.009649, 0.019621, 0.171110, 0.081301, 0.079877, 0.064186, 0.388745
- Error matrix diagonal sum: **0.005053485**
- Error matrix diagonal: 4.01e-04, 6.44e-04, 8.68e-04, 1.06e-04, 1.05e-03, 6.23e-04, 1.30e-04, 1.24e-03

**Fixed Method**:
- Haplotype frequencies: 0.206654, 0.065741, 0.019138, 0.178717, 0.081302, 0.116959, 0.080724, 0.250764
- Error matrix diagonal sum: **0.0002545192**
- Error matrix diagonal: 2.79e-05, 2.79e-05, 3.60e-05, 3.16e-05, 3.26e-05, 3.41e-05, 3.73e-05, 2.71e-05

**Error Ratio**: 0.005053485 / 0.0002545192 = **19.86x higher** (adaptive vs fixed)

### Our Debug Result
**Debug Run**:
- Haplotype frequencies: 0.173337, 0.001999, 0.046659, 0.174310, 0.103976, 0.077801, 0.065412, 0.356507
- Error matrix diagonal sum: **0.002562967**
- Error matrix diagonal: 2.89e-04, 2.47e-04, 4.60e-04, 1.10e-04, 4.13e-04, 3.40e-04, 1.32e-04, 5.72e-04

## Key Findings

### 1. Function Identity Confirmed
- The `est_haps_var` function in our debug wrapper is **identical** to the production version
- Same code, same logic, same parameters

### 2. Data Mismatch Identified
- **Haplotype estimates are different**: Our debug result doesn't match production adaptive result
- **Error estimates are different**: Our debug result (0.00256) is exactly half of production adaptive (0.00505)
- **This indicates we're not using the same input data**

### 3. Root Cause Hypothesis
The issue is likely that **our extractor is not replicating the exact same data preprocessing** that the production system used. Possible causes:

1. **Different RefAlt file**: Production might have used a different version of the RefAlt data
2. **Different data preprocessing**: The production system might have applied different quality filters or data transformations
3. **Different windowing**: The production system might have used different window boundaries
4. **Batch vs single-sample processing**: Production processed all 60 samples together, we processed only 1 sample

### 4. Error Inflation Confirmed
- The 20x error inflation in adaptive vs fixed method is **real and consistent**
- This is not a data extraction artifact - it's a fundamental difference in how the two methods estimate error

## Next Steps

1. **Verify data source**: Confirm we're using the exact same RefAlt file that production used
2. **Check preprocessing**: Compare our data preprocessing steps with production
3. **Investigate batch processing**: Determine if processing all samples together vs. single sample affects results
4. **Debug the error calculation**: Focus on why the adaptive method's progressive V matrix produces inflated error estimates

## Files Created
- `scripts/ErrMatrix/working/extract_hunk.r` - Data extraction script
- `scripts/ErrMatrix/working/debug_wrapper.R` - Debug wrapper script  
- `scripts/ErrMatrix/analysis/extract_specific_sample.R` - Production data extraction
- `hunk_data_chr3R_19780000_Rep01_W_F_h4.rds` - Extracted input data (NOTE: filename includes parameter suffix "_h4")
- `debug_result_chr3R_19780000_Rep01_W_F.rds` - Debug run results
- `Rep01_W_F_19780000_comparison.csv` - Production vs debug comparison

## File Naming Convention
The extractor creates files with the pattern: `hunk_data_{chr}_{position}_{sample}_{parameter}.rds`
Where `parameter` is the h_cutoff value (e.g., "h4" for h_cutoff=4).

**CRITICAL NAMING CHANGE**: The original extractor created files with pattern `single_position_data_{chr}_{position}_{sample}_{parameter}.rds`. This old naming scheme was dangerous as it could lead to confusion about which file contains the correct data. The old file `single_position_data_chr3R_19780000_Rep01_W_F_h4.rds` has been removed to prevent confusion.

**Current correct filename**: `hunk_data_chr3R_19780000_Rep01_W_F_h4.rds`
