# Pair Selection Algorithm - Deep Analysis and Research Summary

## Version History

### v26.0 (Current Baseline): F1 69.81%
- Quality-first triplet selection with adaptive appearance limits
- Multi-phase selection for coverage
- **Stable performance**: 69.81% F1 on dataset 1

### v27.0 (Graph Theory - Greedy Max Coverage): F1 69.28%
- Greedy maximum weighted edge coverage algorithm
- Optimizes for: sum of covered edge weights + diversity bonus
- **Results**:
  - F1: 69.28% (slightly lower than v26)
  - Edge coverage: 328 edges (6% MORE than reference's 309)
  - Total edge weight: 142,574 (17% HIGHER than reference's 122,066)
  - Triplet overlap: 45.9% (same as v26)
  - Dense pair overlap: 72.7% (better than v26)

## Key Discoveries from Analysis

### 1. Reference Selection Pattern (from `find_selection_pattern.py`)
- Reference is **NOT simple top-N selection**
- Only 33.3% of selected triplets are in top-111 quality ranks
- Rank range: 0 to 9,906 out of 31,445 possible triplets
- Average gap between ranks: 90.1
- **Conclusion**: Reference uses **complex global optimization** with multiple criteria

### 2. Edge Coverage vs F1 Score
- v27 achieves **better** edge coverage (328 vs 309) and weight (142k vs 122k)
- BUT F1 is **lower** (69.28% vs 69.81%)
- **Insight**: Maximizing edge coverage alone is insufficient
- Reference algorithm balances multiple objectives beyond just edge weight

### 3. Common vs Different Triplets
From `analyze_missed_triplets.py`:
- Common triplets: avg min_covis=396.3, avg ID span=99.7
- Ref-only triplets: avg min_covis=215.3, avg ID span=545.3
- Gen-only triplets: avg min_covis=386.3, avg ID span=225.5

**Pattern**: Reference includes BOTH high-quality local triplets AND lower-quality bridge triplets

## Why v27 Didn't Improve

### Greedy Algorithm Limitations
1. **Single Objective**: Maximizes edge weight, ignores spatial/geometric constraints
2. **Greedy is Myopic**: Local optimal choices don't guarantee global optimum
3. **Missing Constraints**: Reference may optimize for:
   - Geometric configuration quality
   - Bundle adjustment convergence properties
   - Image graph connectivity patterns
   - Photogrammetric block stability

### Performance Issues
- O(n²) complexity: 536k triplets × 2158 iterations = hours on large datasets
- Need pre-filtering or approximation algorithms

## Hypotheses for Reference Algorithm

### Hypothesis 1: Multi-Objective Optimization
Reference likely uses weighted combination of:
- Edge covisibility quality (weight ~60%)
- Spatial distribution (bridge connections) (weight ~20%)
- Image coverage diversity (weight ~10%)
- Graph connectivity/redundancy (weight ~10%)

### Hypothesis 2: Two-Stage Selection
1. **Stage 1**: Build spanning tree + high-quality local clusters
2. **Stage 2**: Add bridge triplets for global connectivity

### Hypothesis 3: Photogrammetric Constraints
Reference may incorporate domain-specific metrics:
- Convergence angle distribution
- GSD ratio balance
- Optical axis diversity
- Baseline/depth ratios

## Performance Comparison

| Version | F1 Score | Triplet Overlap | Edge Coverage | Edge Weight | Time (DS1) |
|---------|----------|-----------------|---------------|-------------|------------|
| v24     | 69.81%   | 45.9%          | ~256          | ~100k       | <1s        |
| v26     | 69.81%   | 45.9%          | ~256          | ~100k       | <1s        |
| v27     | 69.28%   | 45.9%          | 328           | 142k        | <1s        |

**Conclusion**: Graph theory greedy approach didn't improve F1, though it achieved better coverage metrics.

## Next Approaches to Try

### Approach 1: Reverse Engineering via Machine Learning (Forbidden by user)
- Train classifier on reference selections
- NOT ALLOWED per user request: "不要用机器学习的方法"

### Approach 2: Multi-Objective Optimization
Implement weighted scoring function:
```cpp
score = w1 * edge_quality + w2 * spatial_distribution + w3 * diversity
```
Tune weights using grid search on reference data.

### Approach 3: Hybrid: MST + Quality Clusters + Bridges
1. Build MST for connectivity (guarantees all images connected)
2. Add high-quality local triplets (covis > threshold)
3. Add bridge triplets (large ID span) for long-range connections
4. Balance using appearance constraints

### Approach 4: Graph Partitioning
1. Partition images into spatial regions
2. Select dense triplets within each region
3. Add inter-region bridge triplets
4. Ensures balanced coverage

### Approach 5: Simulated Annealing or Local Search
1. Start with greedy solution
2. Iteratively swap triplets to improve match with reference patterns
3. Use reference statistics as optimization target

## Recommended Next Step

**Try Approach 3: Hybrid MST + Quality Clusters + Bridges**

Rationale:
1. Addresses the "both local AND bridge" pattern we discovered
2. Computationally efficient (no O(n²) loops)
3. Interpretable and tunable
4. Doesn't require machine learning

Implementation plan:
1. Build MST from top-quality edges (ensures connectivity)
2. Find all high-quality complete triangles (min_covis > 200)
3. Select top triangles by quality with appearance limits
4. Add bridge triplets (ID span > 100) if coverage gaps exist
5. Balance to target count (~0.63 × numImages)

Estimated complexity: O(n log n) - fast even for large datasets

## Code Files for Reference

- `test_greedy_coverage.py` - Python prototype of greedy algorithm
- `find_selection_pattern.py` - Analyzes reference selection pattern
- `analyze_missed_triplets.py` - Compares common vs different triplets
- `compare_v26_v27.py` - Detailed v26 vs v27 comparison

## Conclusion

The greedy maximum coverage algorithm (v27) successfully maximized edge coverage and weight but did not improve F1 score, confirming that reference uses multi-objective optimization beyond simple edge quality. The 45.9% triplet overlap plateau suggests we need fundamentally different selection criteria or combinations.

**Current best**: v26.0 with F1 69.81%
**Status**: Need new approach combining structural (MST/bridges) and quality-based selection
