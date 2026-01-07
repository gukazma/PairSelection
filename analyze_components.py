#!/usr/bin/env python3
"""
Analyze reference component structure and image clustering.
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

class UnionFind:
    def __init__(self):
        self.parent = {}
    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return False
        self.parent[px] = py
        return True
    def get_components(self):
        components = defaultdict(set)
        for x in self.parent:
            components[self.find(x)].add(x)
        return list(components.values())

print("Loading data...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])

ref_dense, _, _ = parse_sections('Datas/1/pairs.txt')
ref_set = set(pair_to_tuple(p) for p in ref_dense)

# Build reference components
uf_ref = UnionFind()
for p in ref_set:
    uf_ref.union(p[0], p[1])
ref_components = uf_ref.get_components()

print(f"\nReference has {len(ref_components)} components:")
for i, c in enumerate(sorted(ref_components, key=lambda x: min(x))):
    sorted_ids = sorted(c)
    print(f"\nComponent {i+1}: {len(c)} images")
    print(f"  ID range: {min(c)} - {max(c)}")
    print(f"  Sample IDs: {sorted_ids[:5]}...{sorted_ids[-5:]}")

    # Find cross-component pairs in covisibility
    cross_pairs = []
    for other_comp in ref_components:
        if other_comp == c:
            continue
        for img1 in c:
            for img2 in other_comp:
                pair = (min(img1, img2), max(img1, img2))
                if pair in covis:
                    cross_pairs.append((pair, covis[pair]))

    if cross_pairs:
        cross_pairs.sort(key=lambda x: -x[1])
        print(f"  Best cross-component covis: {cross_pairs[:3]}")

# Analyze why components are separated
print("\n=== Cross-component analysis ===")
all_cross = []
for i, c1 in enumerate(ref_components):
    for j, c2 in enumerate(ref_components):
        if i >= j:
            continue
        for img1 in c1:
            for img2 in c2:
                pair = (min(img1, img2), max(img1, img2))
                if pair in covis:
                    all_cross.append((pair, covis[pair], f"comp{i+1}-comp{j+1}"))

if all_cross:
    all_cross.sort(key=lambda x: -x[1])
    print(f"Total cross-component covisibility pairs: {len(all_cross)}")
    print("Top 10 cross-component pairs (NOT in reference):")
    for p, cv, label in all_cross[:10]:
        in_ref = "IN REF" if p in ref_set else "not in ref"
        print(f"  {p[0]}-{p[1]}: covis={cv}, {label}, {in_ref}")

# Check if these components correspond to ID ranges
print("\n=== Component ID range analysis ===")
# Are components naturally separated by ID?
all_ref_ids = set()
for p in ref_set:
    all_ref_ids.add(p[0])
    all_ref_ids.add(p[1])

sorted_ids = sorted(all_ref_ids)
print(f"All reference IDs: {len(sorted_ids)}")

# Find ID gaps
gaps = []
for i in range(1, len(sorted_ids)):
    gap = sorted_ids[i] - sorted_ids[i-1]
    if gap > 20:
        gaps.append((sorted_ids[i-1], sorted_ids[i], gap))
print(f"\nID gaps > 20:")
for g in gaps[:10]:
    print(f"  {g[0]} -> {g[1]} (gap: {g[2]})")

# Check which MST edges reference DOESN'T have
print("\n=== MST edges analysis ===")
all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])

# Build MST
all_pairs_sorted = sorted(covis.items(), key=lambda x: -x[1])
uf_mst = UnionFind()
mst_edges = []
for pair, cv in all_pairs_sorted:
    if uf_mst.union(pair[0], pair[1]):
        mst_edges.append(pair)
        if len(mst_edges) >= len(all_images) - 1:
            break

mst_set = set(mst_edges)
mst_not_in_ref = mst_set - ref_set
mst_in_ref = mst_set & ref_set

print(f"MST edges: {len(mst_edges)}")
print(f"MST in reference: {len(mst_in_ref)}")
print(f"MST NOT in reference: {len(mst_not_in_ref)}")

print("\nMST edges NOT in reference (would connect components):")
for p in sorted(mst_not_in_ref)[:15]:
    cv = covis.get(p, 0)
    # Check which components this would connect
    comp1 = None
    comp2 = None
    for i, c in enumerate(ref_components):
        if p[0] in c:
            comp1 = i
        if p[1] in c:
            comp2 = i
    print(f"  {p[0]}-{p[1]}: covis={cv}, would connect comp{comp1+1 if comp1 is not None else '?'} to comp{comp2+1 if comp2 is not None else '?'}")
