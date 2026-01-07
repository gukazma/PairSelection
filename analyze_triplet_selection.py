#!/usr/bin/env python3
"""
Analyze why certain images are excluded from reference.
"""

import re
from collections import defaultdict
import glob

def parse_xml_covisibility(xml_path):
    covis = defaultdict(int)
    with open(xml_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    pattern = r'<TiePoint>.*?</TiePoint>'
    for match in re.finditer(pattern, content, re.DOTALL):
        tiepoint = match.group()
        photo_ids = re.findall(r'<PhotoId>(\d+)</PhotoId>', tiepoint)
        photo_ids = [int(x) for x in photo_ids]
        for i in range(len(photo_ids)):
            for j in range(i+1, len(photo_ids)):
                a, b = min(photo_ids[i], photo_ids[j]), max(photo_ids[i], photo_ids[j])
                covis[(a, b)] += 1
    return covis

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

def pair_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

def triplet_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

print("Loading data...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])

ref_dense, _, ref_triplets = parse_sections('Datas/1/pairs.txt')
gen_dense, _, gen_triplets = parse_sections('Datas/1/pairs_generated.txt')

ref_dense_set = set(pair_to_tuple(p) for p in ref_dense)
gen_dense_set = set(pair_to_tuple(p) for p in gen_dense)
ref_tri_set = set(triplet_to_tuple(t) for t in ref_triplets)
gen_tri_set = set(triplet_to_tuple(t) for t in gen_triplets)

# Build adjacency for finding triangles
adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

# Find all possible triangles in covisibility graph
all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])

print(f"\nTotal images: {len(all_images)}")
print(f"Excluded from ref: 2074, 2216, 2271")

# For each excluded image, analyze their covisibility
excluded = [2074, 2216, 2271]
print("\n=== Excluded image analysis ===")
for img in excluded:
    neighbors = adj[img]
    pairs = [(p, covis[p]) for p in covis if img in p]
    pairs.sort(key=lambda x: -x[1])

    print(f"\nImage {img}:")
    print(f"  Neighbors: {len(neighbors)}")
    print(f"  Top covis pairs: {pairs[:5]}")

    # How many triangles does this image participate in?
    triangles = []
    for n1 in neighbors:
        for n2 in neighbors:
            if n1 >= n2:
                continue
            if n2 in adj[n1]:
                tri = tuple(sorted([img, n1, n2]))
                triangles.append(tri)

    print(f"  Participates in {len(triangles)} triangles")

    # Check if any of these triangles are in reference
    in_ref = sum(1 for t in triangles if t in ref_tri_set)
    print(f"  Triangles in reference: {in_ref}")

# Analyze triplets in reference vs generated
print("\n=== Triplet comparison ===")
matching_tri = ref_tri_set & gen_tri_set
only_ref_tri = ref_tri_set - gen_tri_set
only_gen_tri = gen_tri_set - ref_tri_set

print(f"Matching triplets: {len(matching_tri)}")
print(f"Only in ref: {len(only_ref_tri)}")
print(f"Only in gen: {len(only_gen_tri)}")

# What's special about triplets only in reference?
print("\n=== Triplets only in reference ===")
only_ref_list = sorted(only_ref_tri)[:10]
for tri in only_ref_list:
    a, b, c = tri
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    edge_covis = [covis.get(e, 0) for e in edges]
    min_covis = min(edge_covis)
    avg_covis = sum(edge_covis) / 3
    print(f"  {tri}: edge covis = {edge_covis}, min={min_covis}, avg={avg_covis:.0f}")

# What's special about triplets only in generated?
print("\n=== Triplets only in generated ===")
only_gen_list = sorted(only_gen_tri)[:10]
for tri in only_gen_list:
    a, b, c = tri
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    edge_covis = [covis.get(e, 0) for e in edges]
    min_covis = min(edge_covis)
    avg_covis = sum(edge_covis) / 3
    print(f"  {tri}: edge covis = {edge_covis}, min={min_covis}, avg={avg_covis:.0f}")

# What criterion determines if a triplet is in reference?
print("\n=== Triplet selection criterion analysis ===")
# Calculate statistics for ref triplets
ref_stats = []
for tri in ref_tri_set:
    a, b, c = tri
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    edge_covis = [covis.get(e, 0) for e in edges]
    min_covis = min(edge_covis)
    avg_covis = sum(edge_covis) / 3
    id_span = max(tri) - min(tri)
    ref_stats.append((min_covis, avg_covis, id_span))

# Calculate for all possible triplets
all_tri_stats = []
seen = set()
for img in all_images:
    for n1 in adj[img]:
        for n2 in adj[img]:
            if n1 >= n2:
                continue
            if n2 in adj[n1]:
                tri = tuple(sorted([img, n1, n2]))
                if tri in seen:
                    continue
                seen.add(tri)
                a, b, c = tri
                edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
                edge_covis = [covis.get(e, 0) for e in edges]
                min_covis = min(edge_covis)
                avg_covis = sum(edge_covis) / 3
                id_span = max(tri) - min(tri)
                in_ref = tri in ref_tri_set
                all_tri_stats.append((tri, min_covis, avg_covis, id_span, in_ref))

print(f"Total possible triplets: {len(all_tri_stats)}")
print(f"Reference triplets: {len(ref_tri_set)}")

# Sort by various criteria and check reference match
all_tri_stats.sort(key=lambda x: -x[1])  # Sort by min_covis
top_111_min_covis = sum(1 for i, x in enumerate(all_tri_stats[:111]) if x[4])
print(f"\nTop 111 by min edge covis -> in reference: {top_111_min_covis} ({top_111_min_covis*100/111:.1f}%)")

all_tri_stats.sort(key=lambda x: -x[2])  # Sort by avg_covis
top_111_avg_covis = sum(1 for i, x in enumerate(all_tri_stats[:111]) if x[4])
print(f"Top 111 by avg edge covis -> in reference: {top_111_avg_covis} ({top_111_avg_covis*100/111:.1f}%)")

all_tri_stats.sort(key=lambda x: x[3])  # Sort by id_span (smaller is better)
top_111_id_span = sum(1 for i, x in enumerate(all_tri_stats[:111]) if x[4])
print(f"Top 111 by smallest ID span -> in reference: {top_111_id_span} ({top_111_id_span*100/111:.1f}%)")

# Combined score
for w_min, w_avg, w_span in [(1.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.6, 0.3, 0.1), (0.7, 0.2, 0.1)]:
    all_tri_stats.sort(key=lambda x: -(w_min * x[1] + w_avg * x[2] - w_span * x[3] * 0.1))
    top_111_combined = sum(1 for i, x in enumerate(all_tri_stats[:111]) if x[4])
    print(f"Weights (min={w_min}, avg={w_avg}, -span={w_span}): {top_111_combined} ({top_111_combined*100/111:.1f}%)")
