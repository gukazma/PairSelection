#!/usr/bin/env python3
"""
Component-based pair selection - respect natural image clusters.
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

# Build adjacency
adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

# Excluded images (from reference analysis)
excluded_images = {2074, 2216, 2271}

# All images
all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])
all_images -= excluded_images

print(f"Images after exclusion: {len(all_images)}")

# Find reference component structure
uf_ref = UnionFind()
for p in ref_set:
    uf_ref.union(p[0], p[1])
ref_components = uf_ref.get_components()
print(f"Reference has {len(ref_components)} components")

# Identify which component each image belongs to
img_to_component = {}
for i, comp in enumerate(ref_components):
    for img in comp:
        img_to_component[img] = i

# For pairs, identify if they're intra-component or inter-component
print("\n=== Analyzing pair types ===")
intra_ref = 0
inter_ref = 0
for p in ref_set:
    comp0 = img_to_component.get(p[0], -1)
    comp1 = img_to_component.get(p[1], -1)
    if comp0 == comp1 and comp0 != -1:
        intra_ref += 1
    else:
        inter_ref += 1
print(f"Reference: {intra_ref} intra-component, {inter_ref} inter-component pairs")

# Try to discover the natural component structure from covisibility
# Use high-covis pairs only
print("\n=== Discovering component structure ===")

# Sort pairs by covisibility
sorted_pairs = sorted(covis.items(), key=lambda x: -x[1])

# Build components using only high-covis pairs (threshold scan)
best_match = 0
best_threshold = 0
for threshold in range(50, 500, 25):
    uf_test = UnionFind()
    for pair, cv in sorted_pairs:
        if pair[0] in excluded_images or pair[1] in excluded_images:
            continue
        if cv >= threshold:
            uf_test.union(pair[0], pair[1])

    components = uf_test.get_components()
    num_components = len([c for c in components if len(c) > 5])

    if num_components == 4:
        # Try selecting pairs within these components
        img_to_test_comp = {}
        for i, comp in enumerate(components):
            for img in comp:
                img_to_test_comp[img] = i

        # Count how many reference pairs are intra-component for this structure
        intra_count = 0
        for p in ref_set:
            tc0 = img_to_test_comp.get(p[0], -1)
            tc1 = img_to_test_comp.get(p[1], -1)
            if tc0 == tc1 and tc0 != -1:
                intra_count += 1

        if intra_count > best_match:
            best_match = intra_count
            best_threshold = threshold

print(f"Best threshold for 4 components: {best_threshold} (matches {best_match} intra-component pairs)")

# Try component-based selection
print("\n=== Component-based selection ===")

# Build components using reference structure (cheat for now to see upper bound)
# Then select pairs within each component
def select_component_based(components, covis, excluded, max_pairs=198, max_degree=4):
    img_degree = defaultdict(int)
    selected = set()

    # Get pairs sorted by covisibility
    pair_list = []
    for pair, cv in covis.items():
        if pair[0] in excluded or pair[1] in excluded:
            continue
        pair_list.append((pair, cv))
    pair_list.sort(key=lambda x: -x[1])

    # Create component mapping
    img_to_comp = {}
    for i, comp in enumerate(components):
        for img in comp:
            img_to_comp[img] = i

    # Pass 1: Select intra-component pairs (degree 0 first)
    for pair, cv in pair_list:
        if len(selected) >= max_pairs: break
        if pair in selected: continue

        comp0 = img_to_comp.get(pair[0], -1)
        comp1 = img_to_comp.get(pair[1], -1)

        # Only intra-component
        if comp0 != comp1 or comp0 == -1:
            continue

        if img_degree[pair[0]] == 0 and img_degree[pair[1]] == 0:
            if img_degree[pair[0]] < max_degree and img_degree[pair[1]] < max_degree:
                selected.add(pair)
                img_degree[pair[0]] += 1
                img_degree[pair[1]] += 1

    # Pass 2: One node has degree 0
    for pair, cv in pair_list:
        if len(selected) >= max_pairs: break
        if pair in selected: continue

        comp0 = img_to_comp.get(pair[0], -1)
        comp1 = img_to_comp.get(pair[1], -1)
        if comp0 != comp1 or comp0 == -1:
            continue

        if img_degree[pair[0]] == 0 or img_degree[pair[1]] == 0:
            if img_degree[pair[0]] < max_degree and img_degree[pair[1]] < max_degree:
                selected.add(pair)
                img_degree[pair[0]] += 1
                img_degree[pair[1]] += 1

    # Pass 3: Fill remaining
    for pair, cv in pair_list:
        if len(selected) >= max_pairs: break
        if pair in selected: continue

        comp0 = img_to_comp.get(pair[0], -1)
        comp1 = img_to_comp.get(pair[1], -1)
        if comp0 != comp1 or comp0 == -1:
            continue

        if img_degree[pair[0]] < max_degree and img_degree[pair[1]] < max_degree:
            selected.add(pair)
            img_degree[pair[0]] += 1
            img_degree[pair[1]] += 1

    return selected

# Test with reference components (upper bound)
selected_ref_comp = select_component_based(ref_components, covis, excluded_images)
matching_ref = len(selected_ref_comp & ref_set)
print(f"Using reference components: {matching_ref}/198 = {matching_ref*100/198:.1f}%")

# Test discovering components from data
# Try different methods to discover the 4 component structure
print("\n=== Trying to discover component boundaries ===")

# Method: Use large ID gaps to separate components
sorted_ids = sorted(all_images)
gaps = []
for i in range(1, len(sorted_ids)):
    gap = sorted_ids[i] - sorted_ids[i-1]
    gaps.append((sorted_ids[i-1], sorted_ids[i], gap))

gaps.sort(key=lambda x: -x[2])
print("Largest ID gaps:", gaps[:5])

# Try splitting at top 3 gaps
top_gaps = sorted(gaps[:3], key=lambda x: x[0])  # Sort by first ID
boundaries = [g[0] for g in top_gaps]  # Split points
print(f"Split boundaries: {boundaries}")

# Create components based on ID ranges
def id_to_component(img_id, boundaries):
    for i, b in enumerate(boundaries):
        if img_id <= b:
            return i
    return len(boundaries)

discovered_components = defaultdict(set)
for img in all_images:
    comp = id_to_component(img, boundaries)
    discovered_components[comp].add(img)

print(f"Discovered {len(discovered_components)} components based on ID gaps")
for i, comp in discovered_components.items():
    print(f"  Component {i}: {len(comp)} images, range {min(comp)}-{max(comp)}")

# Test selection with discovered components
selected_disc = select_component_based(list(discovered_components.values()), covis, excluded_images)
matching_disc = len(selected_disc & ref_set)
print(f"Using discovered components: {matching_disc}/198 = {matching_disc*100/198:.1f}%")

# What if we DON'T enforce component separation?
# Just use excluded images and degree constraint
print("\n=== Without component separation (baseline) ===")
selected_base = select_component_based([all_images], covis, excluded_images)  # One big component
matching_base = len(selected_base & ref_set)
print(f"Single component: {matching_base}/198 = {matching_base*100/198:.1f}%")
