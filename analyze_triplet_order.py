#!/usr/bin/env python3
"""Analyze reference triplet ordering and patterns"""

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
            if section_idx == 0: dense = section_data
            elif section_idx == 1: refine = section_data
            elif section_idx == 2: triplets = section_data
            section_idx += 1
            i += count + 1
        else:
            i += 1
    return dense, refine, triplets

ref_dense, _, ref_triplets = parse_sections('Datas/1/pairs.txt')

def triplet_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

# Analyze triplet IDs
print("=== Reference triplet ID analysis ===")
triplets = [triplet_to_tuple(t) for t in ref_triplets]
min_ids = [t[0] for t in triplets]
print(f"First triplet min ID: {min_ids[0]}")
print(f"Last triplet min ID: {min_ids[-1]}")
print(f"Min IDs are sorted: {min_ids == sorted(min_ids)}")

# Check if triplets are in order of their minimum ID
print(f"\nFirst 10 triplet min IDs: {min_ids[:10]}")
print(f"Last 10 triplet min IDs: {min_ids[-10:]}")

# Analyze ID ranges/clusters
print("\n=== ID clustering ===")
all_ids = set()
for t in triplets:
    all_ids.update(t)
sorted_ids = sorted(all_ids)
print(f"Total unique IDs: {len(sorted_ids)}")
print(f"ID range: {sorted_ids[0]} to {sorted_ids[-1]}")

# Find gaps in IDs
gaps = []
for i in range(1, len(sorted_ids)):
    gap = sorted_ids[i] - sorted_ids[i-1]
    if gap > 50:
        gaps.append((sorted_ids[i-1], sorted_ids[i], gap))
print(f"\nLarge gaps (>50) in IDs:")
for g in gaps[:5]:
    print(f"  {g[0]} -> {g[1]} (gap: {g[2]})")

# Count triplets per ID range
ranges = [(900, 1300), (1300, 1700), (1700, 2100), (2100, 2500), (2500, 3200)]
print("\n=== Triplets per ID range ===")
for start, end in ranges:
    count = sum(1 for t in triplets if start <= t[0] < end)
    print(f"  {start}-{end}: {count} triplets")
