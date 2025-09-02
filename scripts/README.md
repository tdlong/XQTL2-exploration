# Scripts Directory Organization

This directory contains scripts organized by purpose for the XQTL2 exploration project.

## Directory Structure

### `production/` - Main Production Scripts
These are the scripts you actually use for running analyses:

**Haplotype Estimation:**
- `haplotype_testing_from_table.sh` - Main haplotype testing pipeline
- `haplotype_testing_common.sh` - Common haplotype testing functions
- `haplotype_testing_small_windows.sh` - Small window testing
- `haplotype_estimation_functions.R` - Core haplotype estimation functions
- `run_haplotype_estimation.R` - Haplotype estimation runner

**SNP Imputation:**
- `snp_imputation_from_table.sh` - Main SNP imputation pipeline

**Evaluation & Analysis:**
- `evaluate_imputation_methods.R` - Evaluate different haplotype methods
- `check_snp_imputation_status.R` - Check completion status of imputation jobs
- `create_summary_file_chunked.R` - Create summary files in chunks (preferred)
- `plot_summary_region.R` - Plot specific chromosomal regions
- `summarize_pipeline_results.R` - Summarize overall pipeline results

### `debug/` - Debug and Testing Scripts
One-off scripts for debugging, testing, and development:
- Various test scripts
- Debug scripts
- Comparison scripts
- Development versions
- `create_summary_file.R` - Original summary file creator (superseded by chunked version)

### `archive/` - Old/Legacy Scripts
Outdated or replaced scripts:
- Old REFALT2haps implementations
- Deprecated methods

### `haps2scan/` - Haplotype Scanning
Specialized scripts for haplotype scanning functionality

### `raw2bam2REFALT/` - Data Processing
Scripts for converting raw data to BAM to REFALT format

### `Heterozygosity_tests/` - Heterozygosity Analysis
Scripts for testing heterozygosity patterns

## Usage

For production runs, use scripts from the `scripts/production/` directory (when running from project root):

```bash
# Run haplotype testing
sbatch scripts/production/haplotype_testing_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Run SNP imputation
sbatch scripts/production/snp_imputation_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Check status
Rscript scripts/production/check_snp_imputation_status.R helpfiles/production_slurm_params.tsv process/ZINC2

# Evaluate methods
Rscript scripts/production/evaluate_imputation_methods.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Create summary (preferred)
Rscript scripts/production/create_summary_file_chunked.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Plot specific region
Rscript scripts/production/plot_summary_region.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 870
```

## Notes

- All production scripts now use command-line arguments instead of hard-coded paths
- Scripts are designed to work with any dataset (JUICE, ZINC2, etc.) by specifying parameters
- Debug scripts remain available but are separated from production code
- Archive scripts are kept for reference but should not be used for new analyses
- The `create_summary_file_chunked.R` script is preferred over the original `create_summary_file.R` for memory efficiency
