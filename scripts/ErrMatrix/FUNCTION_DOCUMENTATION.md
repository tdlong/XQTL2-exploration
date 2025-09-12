# Haplotype Error Workbench - Function Documentation

This document describes all functions in `haplotype_error_workbench.R` and their purposes.

## Core Simulation Functions

### `simulate_founders(n_snps, n_founders)`
**Purpose**: Parent-child haplotype simulator generating an n_snps × n_founders (8) 0/1 matrix.

**Concept**:
- Draw K base (unrelated) founders using independent Bernoulli(0.5) SNPs
- The remaining founders are children derived from these bases
- Each child is assigned a target window size W ∈ {300, 750, 1500, 3000}
- For that child, apply per-150-SNP-block flips so that within window W, the child differs from its parent by >20 SNPs on average
- This creates predictable distinguishability timings: smaller W → earlier split

**Parameters**:
- `n_snps`: total number of SNPs (e.g., 3000)
- `n_founders`: number of founders (e.g., 8)

**Returns**: Integer matrix of shape [n_snps, n_founders] with entries in {0,1}

## Haplotype Estimation Functions

### `estimate_haplotypes_list_format_sim(pos, sample_name, df3, founders, h_cutoff, method, window_size_bp, chr, verbose)`
**Purpose**: Core simulated haplotype estimation function with adaptive window algorithm.

**Features**:
- Tests 6 window sizes: 150, 300, 750, 1500, 3000 SNPs
- For each window: filters data, pivots to wide, clusters founders
- Runs LSEI with accumulated constraints from previous windows
- Stops when all founders are distinguishable (8 groups)
- Returns frequency estimates and error matrix
- Uses simulated mode when `pos = -99`

## Wrapper and Simulation Functions

### `run_simulation_wrapper(n_snps, h_cutoff, verbose)`
**Purpose**: Complete simulation wrapper that generates founders + sample → df3, then calls the estimator.

**Process**:
1. Simulates founders with controlled distinguishability
2. Generates sample frequencies using weighted combination + noise
3. Builds df3-like tibble in long format
4. Calls estimator in simulated mode (pos = -99)
5. Prints comparison of true vs estimated haplotype frequencies

## Batch Processing Functions

### `run_batch_df(n_runs, n_snps, h_cutoff)`
**Purpose**: Batch collector returning a tibble with metrics (no printing).

**Returns**: List of results from multiple runs for analysis.

### `run_batch(n_runs, n_snps, h_cutoff)`
**Purpose**: Batch runner with verbose output.

**Features**:
- Runs multiple simulations
- Prints progress and timing information
- Returns list of results

### `run_batch_100_with_summary(h_cutoff)`
**Purpose**: Runs 100 simulations and provides summary statistics.

### `run_batch_df_enhanced(n_runs, n_snps, h_cutoff)`
**Purpose**: Enhanced batch processing with additional metrics and error handling.

### `run_100_with_dataframe(h_cutoff)`
**Purpose**: Runs 100 simulations and returns results as a data frame.

## Analysis and Summary Functions

### `summarize_batch_results(df)`
**Purpose**: Analyzes batch results and computes summary statistics.

**Metrics**:
- Success rates
- Timing statistics
- Group count distributions
- Error matrix properties

### `print_batch_summary(df)`
**Purpose**: Pretty-prints batch results summary.

### `run_one_verbose(run_index, n_snps, h_cutoff, verbose)`
**Purpose**: Runs a single simulation with detailed verbose output.

**Use**: Debugging and detailed analysis of individual runs.

## Data Preparation Functions

### `pre_simulate_data(n_runs, n_snps)`
**Purpose**: Pre-simulates all datasets for benchmarking to avoid timing bias.

**Returns**: List of pre-computed simulation data including:
- Long format data (`df3_long`)
- Wide format data (`df3_wide`)
- Founder matrices
- Sample names and parameters

## Benchmarking Functions

### `benchmark_core_functions(n_runs)`
**Purpose**: Benchmarks the core haplotype estimation function.

**Process**:
1. Pre-simulates all datasets
2. Times the function
3. Calculates performance metrics
4. Estimates billion-scale performance

**Output**:
- Timing results
- Billion-scale time estimates

### `benchmark_lsei_only(n_runs)`
**Purpose**: Micro-benchmark focusing only on LSEI calls to establish theoretical maximum.

**Use**: Identifies if LSEI is the bottleneck and measures pure optimization potential.

## Profiling Functions

### `profile_haplotype_estimator(n_runs)`
**Purpose**: Detailed profiling of haplotype estimation function.

**Features**:
- Line-by-line timing analysis
- Memory usage tracking
- Bottleneck identification
- Optimization recommendations

## Utility Functions

### `pooled_cov(A_full, y_full, groups_vec)`
**Purpose**: Computes pooled covariance at a window for constraint accumulation.

### `fmt_cell_signed(x, diag_cell)`
**Purpose**: Formats covariance matrix cells for compact display.

### `print_V_compact(Vmat)`
**Purpose**: Pretty-prints covariance matrix in compact scientific notation.

## Usage Patterns

### For Development and Testing:
- Use `run_simulation_wrapper()` for single test runs
- Use `run_batch_100_with_summary()` for comprehensive testing
- Use `benchmark_core_functions()` for performance analysis

### For Production Benchmarking:
- Use `pre_simulate_data()` to prepare datasets
- Use `benchmark_core_functions()` for performance analysis
- Use `profile_haplotype_estimator()` for detailed analysis

### For Debugging:
- Use `run_one_verbose()` for detailed single runs
- Use simulated mode (`pos = -99`) to avoid production dependencies

## Key Design Principles

1. **Simulated Mode**: Use `pos = -99` to run in simulation mode without production dependencies
2. **Pre-simulation**: Always pre-simulate data for fair benchmarking
3. **Single Core Function**: Focus on the core `estimate_haplotypes_list_format_sim` function
4. **Result Verification**: Always verify that optimizations don't change results
5. **Scalability**: Estimate performance at billion-scale for production planning