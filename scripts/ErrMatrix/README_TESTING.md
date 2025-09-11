# Wide Format Optimization Testing

## Where to Run From
Run all commands from the **project root directory** (`/Users/tdlong/Desktop/Cursor/XQTL2_exploration/`):

```bash
cd /Users/tdlong/Desktop/Cursor/XQTL2_exploration
```

## Testing Scripts

### 1. Test Production Version
```bash
Rscript scripts/ErrMatrix/test_production.R chr2R process/JUICE helpfiles/JUICE_haplotype_parameters.R
```

### 2. Test Wide Format Version  
```bash
Rscript scripts/ErrMatrix/test_wide_format.R chr2R process/JUICE helpfiles/JUICE_haplotype_parameters.R
```

### 3. Benchmark Both Versions
```bash
Rscript scripts/ErrMatrix/benchmark_wide_vs_production.R chr2R process/JUICE helpfiles/JUICE_haplotype_parameters.R
```

## What Each Script Does

- **test_production.R**: Runs production version in debug mode (100 positions × 1 sample)
- **test_wide_format.R**: Runs wide format version in debug mode (100 positions × 1 sample)  
- **benchmark_wide_vs_production.R**: Runs both versions on same data and compares performance

## Safe Testing
- All scripts use debug mode (limited data)
- Wide format saves to `haplotype_results_list_format_wide_optimized/` (won't overwrite production)
- No risk to existing production data

## Expected Output
- Timing information for each version
- Speedup measurement (hopefully 5-10x faster)
- Verification that results are similar between versions
