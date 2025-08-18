# Adaptive Window Haplotype Estimation Algorithm

## Overview
The adaptive window algorithm estimates founder haplotype frequencies for a given genomic position by progressively expanding the analysis window until founders can be distinguished, then running LSEI estimation at the optimal window size. The algorithm outputs both haplotype frequency estimates and a quality flag indicating reliability.

## Core Algorithm Steps

### 1. Progressive Window Expansion
Start with a small window (e.g., 10kb) and progressively expand to larger windows:
- Window sizes: 10kb → 25kb → 50kb → 100kb → 200kb → 500kb
- Each window is centered on the target position

### 2. Window Size Optimization Loop
For each progressively larger window size:
1. **Extract founder data** for SNPs in the current window
2. **Hierarchical clustering**: Calculate distances and cluster founders  
3. **Distinguishability check**: Cut tree at h_cutoff to count founder groups
4. **Success condition**: If all founders distinguishable (n_groups = n_founders), stop and record optimal window size
5. **Continue expansion**: If not distinguishable, try next larger window size

### 3. LSEI Haplotype Estimation (At Optimal Window)
After finding optimal window size:
1. **Extract data**: Both founder and sample frequencies for SNPs in optimal window
2. **Run LSEI**: Constrained least squares to get B1, B2, ..., AB8 frequencies
3. **Constraints**: Sum to 1, non-negative, lower bound 0.0003
4. **Quality assessment**: Use distinguishability result to set estimate_OK flag

### 4. Single Result Output
For each position/sample combination:
1. **One result only**: Single output per position/sample (not per window tested)
2. **Complete data**: Haplotype frequencies (B1, B2, ..., AB8) from LSEI at optimal window
3. **Quality flag**: estimate_OK (1=reliable, 0=unreliable, NA=LSEI failed)
4. **Metadata**: Final window size, number of SNPs, h_cutoff parameter

## Key Design Principles

### 1. Single Output Per Position/Sample
- **No intermediate results**: Algorithm tests multiple window sizes but outputs only final result
- **Optimal window selection**: Uses smallest window that achieves distinguishability
- **Complete estimation**: Always runs LSEI at chosen window size (if sufficient data)

### 2. Quality Assessment Integration
- **Distinguishability check**: Hierarchical clustering determines if founders can be separated
- **Quality flag**: estimate_OK indicates reliability of haplotype estimates  
- **Flexible threshold**: h_cutoff parameter controls clustering sensitivity

### 3. Efficiency Considerations
- **Early termination**: Stops expanding when distinguishability achieved
- **Data filtering**: Quality filter applied once at data loading stage
- **Memory efficient**: Processes one position/sample at a time

## Mathematical Foundation: Why Constraint Accumulation Works

### The Core Problem: Indistinguishable Founders

When founders are very similar to each other, a fundamental mathematical issue arises:

**LSEI runs perfectly and finds valid solutions, BUT individual founder frequencies become arbitrary.**

#### Example Scenario:
- Founders B1 and B2 cannot be distinguished by clustering
- Their true combined frequency is 0.3
- LSEI will find mathematically valid solutions like:
  - B1=0.1, B2=0.2 (sum = 0.3) ✓
  - B1=0.15, B2=0.15 (sum = 0.3) ✓ 
  - B1=0.05, B2=0.25 (sum = 0.3) ✓

**All solutions are equally valid mathematically, but the individual frequencies are meaningless.**

### Why estimate_OK is Critical

The `estimate_OK` flag captures this fundamental mathematical limitation:

- **`estimate_OK = 1`**: All founders distinguishable → Individual frequencies are meaningful and trustworthy
- **`estimate_OK = 0`**: Some founders indistinguishable → Individual frequencies are arbitrary, only group sums are reliable
- **`estimate_OK = NA`**: LSEI failed → No valid estimates available

### How Constraint Accumulation Solves This

The adaptive window algorithm addresses the indistinguishability problem through progressive constraint accumulation:

#### Small Window (e.g., 10kb):
- B1 and B2 cannot be distinguished
- LSEI estimates: B1=0.1, B2=0.2 (arbitrary individual values)
- **Key insight**: Their sum B1+B2=0.3 is reliable
- **Create constraint**: B1 + B2 = 0.3 (preserve the reliable sum)

#### Larger Window (e.g., 25kb):
- Include accumulated constraint: B1 + B2 = 0.3
- Even if B1 and B2 are still indistinguishable, their sum is constrained correctly
- If other founders become distinguishable, add new constraints
- **Progressive refinement**: Build up reliable group sum constraints

#### Final Window:
- All accumulated constraints ensure group sums remain accurate
- Individual frequencies within indistinguishable groups may still be arbitrary
- But the overall solution respects all reliable group relationships discovered

### The Innovation

**This is why constraint accumulation is the core innovation:**

1. **Preserves reliable information**: Group sums from smaller windows become constraints for larger windows
2. **Prevents information loss**: Reliable relationships discovered early are never lost
3. **Handles partial distinguishability**: Can work even when only some founders are distinguishable
4. **Mathematically sound**: Ensures solutions respect all reliable relationships discovered across window sizes

**Without constraint accumulation, the algorithm would just be a fixed window method with no ability to preserve reliable group relationships discovered at smaller scales.**

## Algorithm Output

### Output Structure
For each genomic position and sample combination:
```R
# Single row with complete information:
chr, pos, sample, h_cutoff, final_window_size, n_snps, estimate_OK, B1, B2, B3, B4, B5, B6, B7, AB8
```

### Output Interpretation
- **Haplotype frequencies**: B1-AB8 values represent the estimated proportion of each founder haplotype
- **estimate_OK values**:
  - `1`: Reliable estimates (founders distinguishable at h_cutoff)
  - `0`: Unreliable estimates (founders not distinguishable, but LSEI ran)
  - `NA`: No estimates available (LSEI failed due to insufficient data)
- **final_window_size**: Optimal window size used for estimation (in bp)
- **n_snps**: Number of SNPs used in the final estimation

## Implementation Notes

- **Single result per position/sample**: No intermediate results from window testing phases
- **Optimal window selection**: Uses smallest window that achieves distinguishability
- **Complete data requirement**: Both founder and sample data must be present
- **Quality filtering**: Applied once at data loading stage for efficiency
- **Dynamic founder count**: Works with any number of founders defined in parameter file

### 6. Convergence Criteria
The algorithm stops when:
- **All founders are separated** (n_groups == n_founders) - **OPTIMAL OUTCOME**
- **Maximum window size reached** (500kb) without full separation - **SUBOPTIMAL BUT VALID**
- **No meaningful group changes** across multiple window sizes

## Key Parameters

### h_cutoff
- Controls the "height" at which the hierarchical clustering tree is cut
- Lower values (e.g., 4) = more aggressive grouping = fewer founder groups
- Higher values (e.g., 10) = less aggressive grouping = more founder groups
- Should produce different results for different values

### Window Sizes
- Start small (10kb) for better founder separation
- Expand progressively to capture more SNPs
- Balance between statistical power and founder distinguishability

## Expected Behavior

### With Different h_cutoff Values
- **h4**: More aggressive grouping → fewer founder groups → more constrained estimates
- **h10**: Less aggressive grouping → more founder groups → less constrained estimates
- **Result**: Different founder frequency estimates for different h_cutoff values

### Progressive Improvement
- **Smaller windows**: Better founder separation, less statistical power
- **Larger windows**: More statistical power, potentially different tree structure
- **Tree adaptation**: Structure may change (reshuffling) as more SNPs are included
- **Accumulated constraints**: Preserve good separation from smaller windows

## Algorithm Advantages

1. **Adaptive**: Window size adapts to local genomic structure
2. **Constraint-preserving**: Good founder separation from small windows is preserved
3. **Statistically robust**: Larger windows provide more data for estimation
4. **Parameter-controlled**: h_cutoff allows tuning of grouping aggressiveness
5. **Resilient**: Handles tree reshuffling gracefully as statistical power increases

## Implementation Notes

### Critical Components
1. **Hierarchical clustering** must be performed at each window
2. **Smart constraint detection** must identify meaningful group changes
3. **Constraint accumulation** must preserve group and individual constraints
4. **LSEI solver** must handle the accumulated constraint matrix
5. **Progress tracking** must monitor group composition across windows

### Common Pitfalls
1. **Redundant constraints**: Adding sum=1 constraint when all founders are grouped
2. **Unnecessary LSEI calls**: Running solver when groups haven't changed
3. **Wrong constraint formulation**: Not properly translating groups to constraints
4. **Ignoring reshuffling**: Not adapting to tree structure changes

### Handling Edge Cases
1. **No separation achieved**: Return best result from largest window
2. **Tree reshuffling**: Adapt constraints and continue
3. **Insufficient SNPs**: Skip window or return NA
4. **LSEI failures**: Continue to next window size

## Current Status

### Bug Identified and Fixed ✅
**The redundant constraint bug has been identified and fixed in the production code:**

1. **Root Cause**: The 10kb window (1 founder group) was adding a redundant "sum=1" constraint that was already enforced by the base constraint system.

2. **The Problem**: This redundant constraint caused matrix singularity and LSEI failures, leading to identical results across all h_cutoff values.

3. **The Fix Implemented**: 
   - Added `groups_changed()` function to detect meaningful group composition changes
   - Only run LSEI and accumulate constraints when groups meaningfully change
   - Skip redundant constraint accumulation when all founders are in one group

4. **Verification**: Test script confirmed the fix works - h4 vs h10 now produce different results with frequency differences of 0.055-0.118.

### Production Code Status
**The current production version (`scripts/REFALT2haps.AdaptWindow.Single.R`) now has the fixed constraint logic:**

- ✅ **Full hierarchical clustering algorithm** implemented
- ✅ **Progressive window expansion** working correctly
- ✅ **Smart constraint accumulation** - only when groups meaningfully change
- ✅ **Redundant constraint bug** fixed
- 🔄 **Pipeline re-execution** in progress with fixed code

### Expected Results
With the bug fix, the algorithm should now:
- **Produce different results** for different h_cutoff values (h4, h6, h8, h10)
- **Show meaningful clustering** effects across the parameter range
- **Process full chromosome** coverage (not just 2429 positions)
- **Demonstrate adaptive behavior** as intended

### Next Steps
1. **Monitor current pipeline** - Adaptive window jobs are running with fixed code
2. **Verify different performance** across h_cutoff values once pipeline completes
3. **Extend to other chromosomes** using the working algorithm
4. **Parameter optimization** - find optimal h_cutoff values for JUICE dataset
