#!/usr/bin/env python3
"""
Analyze: What if we select pairs first, then form triplets?

Reference algorithm might work as:
1. Select ~198 high-quality pairs (with image coverage constraint)
2. Find triplets that are fully covered by those pairs
"""

from collections import defaultdict

def parse_sections(filepath):
    dense, refine, triplets = [], [], []
    with open(filepath, 'r') as f:
        lines = [l.strip() for l in f.readlines()]

    i, section_idx = 0, 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('#') or line == '':
            i += 1
            continue
        if line.isdigit():
            count = int(line)
            section_data = []
            for j in range(i+1, min(i+1+count, len(lines))):
                if lines[j]:
                    section_data.append(lines[j])
            if section_idx == 0:
                dense = section_data
            elif section_idx == 1:
                refine = section_data
            elif section_idx == 2:
                triplets = section_data
            section_idx += 1
            i += count + 1
        else:
            i += 1
    return dense, refine, triplets

# Parse reference
ref_dense, ref_refine, ref_triplets = parse_sections('Datas/1/pairs.txt')

print(f"Reference: {len(ref_dense)} dense, {len(ref_triplets)} triplets")

# Convert to sets
def pair_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

def triplet_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

dense_set = set(pair_to_tuple(p) for p in ref_dense)
triplet_set = set(triplet_to_tuple(t) for t in ref_triplets)

# Check: are all triplet edges in dense pairs?
all_triplet_edges = set()
for t in triplet_set:
    a, b, c = t
    all_triplet_edges.add((min(a,b), max(a,b)))
    all_triplet_edges.add((min(a,c), max(a,c)))
    all_triplet_edges.add((min(b,c), max(b,c)))

print(f"All triplet edges: {len(all_triplet_edges)}")
print(f"Triplet edges in dense: {len(all_triplet_edges & dense_set)}")
print(f"Dense pairs not in triplets: {len(dense_set - all_triplet_edges)}")

# So we know all dense pairs are triplet edges
# Now: from dense pairs, what triplets can we form?
# A triplet exists if all 3 edges are in dense_set

# Build adjacency from dense pairs
adj = defaultdict(set)
for p in dense_set:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

# Find all possible triplets from dense pairs
possible_triplets = set()
for edge in dense_set:
    a, b = edge
    common = adj[a] & adj[b]
    for c in common:
        tri = tuple(sorted([a, b, c]))
        # Check all 3 edges exist
        if (min(a,b), max(a,b)) in dense_set and \
           (min(a,c), max(a,c)) in dense_set and \
           (min(b,c), max(b,c)) in dense_set:
            possible_triplets.add(tri)

print(f"\nTriplets that can be formed from dense pairs: {len(possible_triplets)}")
print(f"Reference triplets: {len(triplet_set)}")
print(f"Matching: {len(possible_triplets & triplet_set)}")

# So reference triplets are a SUBSET of possible triplets from dense pairs
# This confirms: pairs first, then select triplets

# Now let's check which triplets are selected and why
# Reference has 111 triplets but we can form more from 198 pairs

# Analyze: each reference triplet - how many of its edges are in dense?
for t in list(triplet_set)[:5]:
    a, b, c = t
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    in_dense = sum(1 for e in edges if e in dense_set)
    print(f"Triplet {t}: {in_dense}/3 edges in dense")
