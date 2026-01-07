#!/usr/bin/env python3
"""
Analyze network structure of reference pairs to understand bridge pair selection.
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

# Union-Find for connected components
class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return False
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1
        return True

    def get_components(self):
        components = defaultdict(set)
        for x in self.parent:
            components[self.find(x)].add(x)
        return list(components.values())

print("Loading data...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])
print(f"Loaded {len(covis)} covisibility pairs")

ref_dense, _, ref_triplets = parse_sections('Datas/1/pairs.txt')
gen_dense, _, gen_triplets = parse_sections('Datas/1/pairs_generated.txt')

ref_set = set(pair_to_tuple(p) for p in ref_dense)
gen_set = set(pair_to_tuple(p) for p in gen_dense)

# Analyze connected components
print("\n=== Connected component analysis ===")

# All images with covisibility
all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])
print(f"Total images with covisibility: {len(all_images)}")

# Check if reference pairs form a connected graph
uf_ref = UnionFind()
for p in ref_set:
    uf_ref.union(p[0], p[1])
ref_components = uf_ref.get_components()
print(f"Reference pairs form {len(ref_components)} connected component(s)")
if len(ref_components) <= 5:
    for i, c in enumerate(sorted(ref_components, key=len, reverse=True)[:5]):
        print(f"  Component {i+1}: {len(c)} images, ID range: {min(c)}-{max(c)}")

# Check generated
uf_gen = UnionFind()
for p in gen_set:
    uf_gen.union(p[0], p[1])
gen_components = uf_gen.get_components()
print(f"Generated pairs form {len(gen_components)} connected component(s)")

# Identify bridge pairs in reference (pairs that connect different ID regions)
print("\n=== Bridge pair analysis ===")

# Define ID regions based on gaps
sorted_ids = sorted(all_images)
gaps = []
for i in range(1, len(sorted_ids)):
    gap = sorted_ids[i] - sorted_ids[i-1]
    if gap > 100:
        gaps.append((sorted_ids[i-1], sorted_ids[i], gap))

print(f"Large gaps (>100) in image IDs:")
for g in gaps:
    print(f"  {g[0]} -> {g[1]} (gap: {g[2]})")

# Find pairs that cross these gaps
def crosses_gap(pair, gaps):
    a, b = pair
    for g in gaps:
        if a <= g[0] and b >= g[1]:
            return True
    return False

ref_bridge = [p for p in ref_set if crosses_gap(p, gaps)]
gen_bridge = [p for p in gen_set if crosses_gap(p, gaps)]
print(f"\nPairs crossing major gaps - Ref: {len(ref_bridge)}, Gen: {len(gen_bridge)}")

if ref_bridge:
    print("Reference bridge pairs:")
    for p in sorted(ref_bridge)[:10]:
        print(f"  {p[0]}-{p[1]}, covis={covis.get(p, 0)}, span={p[1]-p[0]}")

# Analyze what makes reference bridge pairs special
print("\n=== Reference bridge pair selection criteria ===")
# For each bridge pair, what other pairs exist for those images?
for bp in sorted(ref_bridge)[:5]:
    img1, img2 = bp
    print(f"\nBridge pair {img1}-{img2} (covis={covis.get(bp, 0)}):")

    # Find all pairs involving img1 and img2
    img1_pairs = [(p, covis[p]) for p in covis if img1 in p and p != bp]
    img2_pairs = [(p, covis[p]) for p in covis if img2 in p and p != bp]

    img1_pairs.sort(key=lambda x: -x[1])
    img2_pairs.sort(key=lambda x: -x[1])

    print(f"  Image {img1} has {len(img1_pairs)} other pairs, top covis: {[x[1] for x in img1_pairs[:3]]}")
    print(f"  Image {img2} has {len(img2_pairs)} other pairs, top covis: {[x[1] for x in img2_pairs[:3]]}")

    # Is this bridge pair the only connection between regions?
    img1_region = set(p[0] if p[1] == img1 else p[1] for p, _ in img1_pairs)
    img2_region = set(p[0] if p[1] == img2 else p[1] for p, _ in img2_pairs)
    overlap = img1_region & img2_region
    print(f"  Overlap between neighbor sets: {len(overlap)}")

# Check minimum spanning tree properties
print("\n=== Minimum Spanning Tree analysis ===")
# Sort all covisibility pairs by covis (descending) and find MST
all_pairs_sorted = sorted(covis.items(), key=lambda x: -x[1])
uf_mst = UnionFind()
mst_pairs = []
for pair, cv in all_pairs_sorted:
    if uf_mst.union(pair[0], pair[1]):
        mst_pairs.append(pair)
        if len(mst_pairs) >= len(all_images) - 1:
            break

mst_set = set(mst_pairs)
mst_in_ref = len(mst_set & ref_set)
mst_in_gen = len(mst_set & gen_set)
print(f"MST edges: {len(mst_pairs)}")
print(f"MST edges in reference: {mst_in_ref} ({mst_in_ref*100/len(mst_pairs):.1f}%)")
print(f"MST edges in generated: {mst_in_gen} ({mst_in_gen*100/len(mst_pairs):.1f}%)")

# Which MST edges are missing from reference?
mst_not_in_ref = mst_set - ref_set
print(f"\nMST edges NOT in reference: {len(mst_not_in_ref)}")
for p in sorted(mst_not_in_ref)[:5]:
    print(f"  {p[0]}-{p[1]}, covis={covis.get(p, 0)}")

# Check if reference has all necessary connectivity
print("\n=== Connectivity verification ===")
uf_check = UnionFind()
for p in ref_set:
    uf_check.union(p[0], p[1])
covered_images = set()
for p in ref_set:
    covered_images.add(p[0])
    covered_images.add(p[1])
uncovered = all_images - covered_images
print(f"Images covered by reference: {len(covered_images)}/{len(all_images)}")
if uncovered:
    print(f"Uncovered images: {sorted(uncovered)[:10]}...")
