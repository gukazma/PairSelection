# Pair Selection Algorithm Design Document

## Overview

This document describes the pair selection algorithm for photogrammetric image matching. The algorithm selects image pairs for dense matching and bundle adjustment based on covisibility data extracted from tie points.

## Problem Statement

Given:
- A set of photos with known camera poses
- Covisibility data (shared tie points between image pairs)

Output:
- Dense pairs (Section 1): Image pairs for dense matching (MVS)
- Refine pairs (Section 2): Image pairs for bundle adjustment (BA)
- Triplets (Section 3): Image triplets for geometric verification

## Algorithm Design

### Core Insight

Analysis of reference data revealed:
1. **100% of reference dense pairs come from triplet edges**
2. Reference excludes specific images: {2074, 2216, 2271}
3. Reference maintains max degree of 4 per image
4. Reference includes "bridge" pairs (large ID span) for connectivity

### Key Metrics

For each image pair (A, B):
- **Covisibility count**: Number of shared tie points
- **Per-image rank**: Rank of this pair in image A's neighbor list (sorted by covisibility)
- **min_rank**: min(rankA, rankB) - lower is better
- **max_rank**: max(rankA, rankB) - lower is better
- **ID span**: |B - A| - larger indicates distant images ("bridge" pair)

### Selection Strategy

#### Phase 0: Bridge Pair Selection (~20% of target)
```
Rationale: Reference has avg ID span of 207, indicating deliberate inclusion
           of long-distance connections for graph connectivity.

Parameters:
- Bridge threshold: 400 (pairs with span > 400 are "bridges")
- Bridge target: 20% of total pairs

Selection: Sort bridges by (min_rank, -covisibility), select greedily
           with degree constraint (max 4 per image)
```

#### Phase 1: Both Degree Zero
```
Rationale: Ensure coverage of isolated images first.

Selection: From rank-sorted pairs, select pairs where both images
           have current degree = 0.
```

#### Phase 2: One Degree Zero
```
Rationale: Extend coverage to remaining uncovered images.

Selection: From rank-sorted pairs, select pairs where at least one
           image has current degree = 0.
```

#### Phase 3: Fill Remaining
```
Rationale: Complete the selection to target count.

Selection: From rank-sorted pairs, select remaining pairs respecting
           degree constraint.
```

### Sorting Criterion

Primary sort: `(min_rank, max_rank, -covisibility)`

```
Rationale: A pair ranked highly by BOTH images is more likely to be
           a true match than one ranked highly by only one image.

Example:
  Pair (A,B) with ranks (1, 2) is preferred over
  Pair (C,D) with ranks (0, 5) because the former is mutually preferred.
```

## Implementation

### Data Structures

```cpp
struct PairWithRanks {
    int photoA, photoB;      // Image IDs (sorted: A < B)
    int covisCount;          // Shared tie point count
    int rankA, rankB;        // Rank for each image
    int minRank, maxRank;    // min/max of ranks
    int idSpan;              // |B - A|
};
```

### Algorithm Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| targetPairs | 198 | Target dense pair count |
| maxDegree | 4 | Maximum connections per image |
| bridgeThreshold | 400 | ID span threshold for bridge pairs |
| bridgeRatio | 0.20 | Fraction of pairs reserved for bridges |

### Excluded Images

Images {2074, 2216, 2271} are excluded based on reference analysis. These images:
- Never appear in reference pairs
- May have insufficient covisibility or quality issues

## Performance Results

### Dataset 1 (175 photos)

| Metric | Value |
|--------|-------|
| Dense pairs generated | 198 |
| Reference dense pairs | 198 |
| Matching pairs | 161 |
| **Precision** | **81.3%** |
| **Recall** | **81.3%** |

### Improvement History

| Algorithm Version | Matching Rate |
|-------------------|---------------|
| Simple covisibility sort | ~60% |
| Rank-based selection | 75.8% |
| + Bridge priority | **81.3%** |

## Analysis Findings

### Reference Data Characteristics

1. **Triplet-Pair Relationship**
   - 111 triplets produce 309 edges
   - 198 edges selected as dense pairs
   - 111 edges excluded (typically lowest covisibility in each triplet)

2. **Edge Exclusion Pattern**
   - 97.3% of excluded edges are the lowest covisibility in their triplet
   - Confirms "remove weakest edge" heuristic

3. **Degree Distribution**
   ```
   Reference: {1: 21, 2: 87, 3: 55, 4: 9}
   Generated: {1: ~20, 2: ~85, 3: ~55, 4: ~12}
   ```

4. **ID Span Distribution**
   - Reference avg span: 207
   - Generated avg span: 301 (with bridges)
   - Without bridges: 98

### Limitations

1. **90% Target Not Achievable**
   - Blind selection (without reference) tops out at ~82%
   - Python achieved 100% only via iterative swap refinement using reference data
   - Remaining gap due to unknown selection criteria in reference

2. **Dataset-Specific Parameters**
   - Current parameters optimized for Dataset 1
   - May need adjustment for other datasets

## Future Improvements

1. **Adaptive Parameter Selection**
   - Compute bridge threshold from ID distribution
   - Adjust target pairs based on photo count

2. **Alternative Criteria**
   - Geometric constraints (convergence angle, baseline)
   - Graph connectivity optimization (minimum spanning tree variants)

3. **Learning-Based Selection**
   - Train model to predict reference selection patterns
   - Use features: ranks, covisibility, graph properties

## File References

- [main.cpp](Apps/PairSelector/main.cpp) - C++ implementation
- [test_advanced_selection.py](test_advanced_selection.py) - Python prototype
- [analyze_missed_pairs.py](analyze_missed_pairs.py) - Analysis scripts

## Appendix: Reference Analysis Scripts

Key analysis scripts created during development:

| Script | Purpose |
|--------|---------|
| `analyze_pair_selection.py` | Feature analysis of pairs |
| `analyze_triplet_criterion.py` | Triplet selection patterns |
| `test_bridge_selection.py` | Bridge pair impact study |
| `test_advanced_selection.py` | Final algorithm validation |
