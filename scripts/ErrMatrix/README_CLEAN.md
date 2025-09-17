# ErrMatrix Directory - Clean Structure

## Working Scripts (`working/`)
- `extract_hunk.r` - Extract single position data (self-contained)
- `run_estimator_locally.R` - Run estimator on extracted data
- `compare_hunk_vs_production.R` - Compare hunk results vs production results

## Analysis Scripts (`analysis/`)
- `compare_hap_and_err_diffs.R` - Compare adaptive vs fixed methods
- `analyze_testing_positions.R` - Analyze testing positions data

## Core Files
- `BASE_VAR_WIDE.R` - Production haplotype estimator
- `haplotype_error_workbench.R` - Error analysis tools
- `run_all_chroms.slurm` - SLURM job script

## Archive (`archive/`)
- Historical versions and experimental code

## Workflow
1. Run `extract_hunk.r` on cluster to get single position data
2. Download `hunk_data_*.rds` file locally
3. Run `compare_hunk_vs_production.R` to compare results
4. Use analysis scripts to investigate differences
