# Adaptive Window Haplotype Estimation Algorithm

## Overview
The adaptive window algorithm is designed to estimate founder haplotype frequencies for a given genomic position by progressively expanding the analysis window while accumulating constraints from smaller windows. The key insight is that smaller windows may provide better founder separation, and these constraints can be carried forward to larger windows.

## Core Algorithm Steps

### 1. Progressive Window Expansion
Start with a small window (e.g., 10kb) and progressively expand to larger windows:
- Window sizes: 10kb → 25kb → 50kb → 100kb → 200kb → 500kb
- Each window is centered on the target position

### 2. Hierarchical Clustering at Each Window
For each window size:
1. **Extract founder genotype data** for all SNPs in the window
2. **Calculate pairwise distances** between founders using `dist(t(founder_matrix))`
3. **Perform hierarchical clustering** using `hclust()` 
4. **Cut the tree** at the specified `h_cutoff` using `cutree(hclust_result, h = h_cutoff)`
5. **Identify founder groups** based on the clustering

### 3. Smart Constraint Accumulation
The key innovation is carrying forward constraints only when meaningful changes occur:
- **Detect meaningful group changes**: Not just count, but check if group composition changed
- **Skip LSEI when groups unchanged**: Reuse previous constraints and results
- **Accumulate constraints progressively**: As groups split or reshuffle
- **Handle reshuffling gracefully**: Tree structure may change with more SNPs

### 4. Constraint Types
Two types of constraints are accumulated:

#### Group Constraints
For founder groups that cannot be separated in larger windows:
```
Group 1: B1 + B2 + B3 = 0.45  (from smaller window)
Group 2: B4 + B5 = 0.30       (from smaller window)
```

#### Individual Constraints  
For single founders that were isolated:
```
B6 = 0.15  (from smaller window)
B7 = 0.08  (from smaller window)
B8 = 0.02  (from smaller window)
```

### 5. Constrained Least Squares Estimation
At each window, solve:
```
minimize ||A*x - b||^2
subject to:
  E*x = F  (equality constraints from accumulated constraints)
  G*x >= H (inequality constraints: all frequencies >= 0.0003)
  sum(x) = 1 (frequencies sum to 1)
```

Where:
- `A` = founder genotype matrix for current window
- `b` = sample frequencies for current window
- `E, F` = accumulated constraints from smaller windows
- `G, H` = non-negativity constraints

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

### Code Analysis Results
After examining the git history and current code:

1. **The current production version (`scripts/REFALT2haps.AdaptWindow.Single.R`) DOES implement the full hierarchical clustering algorithm** with:
   - Progressive window expansion
   - Hierarchical clustering using `hclust()` and `cutree()`
   - Constraint accumulation from smaller windows
   - Proper LSEI solving with accumulated constraints

2. **The August 12th version (`scripts/REFALT2haps.AdaptWindow.R`) was the original chromosome-wide scanner** that tested multiple h_cutoff values in a single run.

3. **The current version is the single-parameter version** that runs one h_cutoff value at a time, but uses the same core algorithm.

### The Mystery
Despite having the correct algorithm, the cluster results show:
- Identical founder frequencies across all h_cutoff values (h4, h6, h8, h10)
- Always 8 founder groups (no clustering effect)
- Only 2429 positions processed (subset, not full chromosome)

This suggests the issue might be:
1. **Data-specific**: The founder genotypes might be too similar for h_cutoff to matter
2. **Parameter range**: The h_cutoff values (4, 6, 8, 10) might be too similar
3. **Implementation bug**: Something subtle in the constraint accumulation
4. **Different version run**: The cluster might have run a different version

### Recent Debugging Insights
The verbose test script revealed:
1. **The algorithm IS working correctly** - it progresses from 1 group (10kb) to 7 groups (25kb) to 8 groups (50kb)
2. **A bug exists**: The 10kb window (1 group) adds a redundant constraint (sum=1) which is already enforced
3. **This redundant constraint** could cause matrix singularity and LSEI failures
4. **The fix**: Only accumulate constraints when groups meaningfully change, not just when they exist

## Testing Strategy
The `scripts/test_working_adaptive.R` script tests this algorithm on a couple of positions with different h_cutoff values to verify that:
1. Hierarchical clustering is working
2. Different h_cutoff values produce different results
3. Constraint accumulation is functioning properly
4. The algorithm converges appropriately
5. Redundant constraints are not accumulated

This will help determine if the issue is with the algorithm itself or with the specific data/parameters used on the cluster.
