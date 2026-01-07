#!/usr/bin/env python3
"""
Deep analysis of triplet selection criterion in reference data.
Goal: Understand exactly how reference selects 111 triplets from all possible triangles.
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

print("Loading data...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])

ref_dense, _, ref_triplets = parse_sections('Datas/1/pairs.txt')

# Build adjacency
adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

# All images in covisibility
all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])

# Excluded images
excluded_images = {2074, 2216, 2271}
valid_images = all_images - excluded_images

print(f"Total images: {len(all_images)}")
print(f"Valid images (after exclusion): {len(valid_images)}")

# Parse reference triplets
ref_tri_set = set()
for t in ref_triplets:
    ids = tuple(sorted(map(int, t.split())))
    ref_tri_set.add(ids)

print(f"\nReference triplets: {len(ref_tri_set)}")

# Find ALL possible triplets (triangles in covisibility graph)
print("\nFinding all possible triplets...")
all_triplets = []
seen = set()
for img in valid_images:
    neighbors = [n for n in adj[img] if n in valid_images]
    for i, n1 in enumerate(neighbors):
        for n2 in neighbors[i+1:]:
            if n2 in adj[n1]:
                tri = tuple(sorted([img, n1, n2]))
                if tri not in seen:
                    seen.add(tri)
                    # Calculate triplet features
                    a, b, c = tri
                    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
                    edge_covis = [covis.get(e, 0) for e in edges]
                    min_covis = min(edge_covis)
                    max_covis = max(edge_covis)
                    sum_covis = sum(edge_covis)
                    avg_covis = sum_covis / 3
                    id_span = max(tri) - min(tri)
                    in_ref = tri in ref_tri_set

                    all_triplets.append({
                        'tri': tri,
                        'edges': edges,
                        'edge_covis': edge_covis,
                        'min_covis': min_covis,
                        'max_covis': max_covis,
                        'sum_covis': sum_covis,
                        'avg_covis': avg_covis,
                        'id_span': id_span,
                        'in_ref': in_ref
                    })

print(f"Total possible triplets: {len(all_triplets)}")
ref_count = sum(1 for t in all_triplets if t['in_ref'])
print(f"Reference triplets found: {ref_count}")

# Analyze differences between ref and non-ref triplets
ref_tris = [t for t in all_triplets if t['in_ref']]
non_ref_tris = [t for t in all_triplets if not t['in_ref']]

def stats(tris, key):
    vals = [t[key] for t in tris]
    return min(vals), max(vals), sum(vals)/len(vals) if vals else 0

print("\n=== Feature Comparison ===")
for feat in ['min_covis', 'max_covis', 'sum_covis', 'avg_covis', 'id_span']:
    ref_s = stats(ref_tris, feat)
    non_s = stats(non_ref_tris, feat)
    print(f"\n{feat}:")
    print(f"  Reference: min={ref_s[0]:.1f}, max={ref_s[1]:.1f}, avg={ref_s[2]:.1f}")
    print(f"  Non-ref:   min={non_s[0]:.1f}, max={non_s[1]:.1f}, avg={non_s[2]:.1f}")

# Try different selection criteria
print("\n=== Testing Selection Criteria ===")

def test_criterion(triplets, score_fn, name, target=111):
    sorted_tris = sorted(triplets, key=score_fn, reverse=True)
    selected = sorted_tris[:target]
    matching = sum(1 for t in selected if t['in_ref'])
    print(f"{name}: {matching}/{target} = {matching*100/target:.1f}%")
    return selected

# Simple criteria
test_criterion(all_triplets, lambda t: t['min_covis'], "min_edge_covis")
test_criterion(all_triplets, lambda t: t['sum_covis'], "sum_covis")
test_criterion(all_triplets, lambda t: t['avg_covis'], "avg_covis")
test_criterion(all_triplets, lambda t: -t['id_span'], "smallest_id_span")

# Combined criteria
test_criterion(all_triplets, lambda t: t['min_covis'] * 2 + t['avg_covis'], "2*min + avg")
test_criterion(all_triplets, lambda t: t['min_covis'] + t['sum_covis']/100, "min + sum/100")

# Image degree constraint
print("\n=== With Image Coverage Constraint ===")

def select_with_coverage(triplets, score_fn, target=111, max_per_image=3):
    """Select triplets ensuring each image appears in limited triplets"""
    sorted_tris = sorted(triplets, key=score_fn, reverse=True)
    selected = []
    img_count = defaultdict(int)

    for t in sorted_tris:
        if len(selected) >= target:
            break
        tri = t['tri']
        # Check if any image has too many triplets
        if all(img_count[img] < max_per_image for img in tri):
            selected.append(t)
            for img in tri:
                img_count[img] += 1

    matching = sum(1 for t in selected if t['in_ref'])
    print(f"Score: {score_fn.__name__ if hasattr(score_fn, '__name__') else 'custom'}, max_per_img={max_per_image}: {matching}/{len(selected)} = {matching*100/len(selected):.1f}%")
    return selected

for max_img in [2, 3, 4, 5]:
    select_with_coverage(all_triplets, lambda t: t['min_covis'], max_per_image=max_img)

# Analyze images in reference triplets
print("\n=== Image Triplet Participation ===")
ref_img_count = defaultdict(int)
for t in ref_tris:
    for img in t['tri']:
        ref_img_count[img] += 1

counts = list(ref_img_count.values())
print(f"Images in reference triplets: {len(ref_img_count)}")
print(f"Triplet count per image: min={min(counts)}, max={max(counts)}, avg={sum(counts)/len(counts):.1f}")
count_dist = defaultdict(int)
for c in counts:
    count_dist[c] += 1
print(f"Distribution: {dict(sorted(count_dist.items()))}")

# Check if reference uses specific images
print("\n=== Looking for Image Selection Pattern ===")
# Images in ref triplets vs all possible triplets
all_tri_images = set()
for t in all_triplets:
    for img in t['tri']:
        all_tri_images.add(img)

ref_tri_images = set(ref_img_count.keys())
print(f"Images in all possible triplets: {len(all_tri_images)}")
print(f"Images in reference triplets: {len(ref_tri_images)}")
print(f"Images only in ref: {len(ref_tri_images - all_tri_images)}")
print(f"Images not in ref triplets: {len(all_tri_images - ref_tri_images)}")

not_in_ref = all_tri_images - ref_tri_images
if not_in_ref:
    print(f"  Excluded from ref triplets: {sorted(not_in_ref)[:20]}")

# Check reference pair coverage
print("\n=== Reference Dense Pair Analysis ===")
ref_dense_set = set()
for p in ref_dense:
    ids = tuple(sorted(map(int, p.split())))
    ref_dense_set.add(ids)

# From selected triplets, extract pairs
def extract_pairs_from_triplets(triplets):
    """Extract 2 best edges from each triplet (remove lowest covis)"""
    pairs = set()
    for t in triplets:
        edges = t['edges']
        edge_covis = t['edge_covis']
        # Sort edges by covis, remove lowest
        sorted_edges = sorted(zip(edges, edge_covis), key=lambda x: -x[1])
        # Keep top 2
        for edge, _ in sorted_edges[:2]:
            pairs.add(edge)
    return pairs

# Test triplet-first approach
print("\n=== Triplet-First Algorithm Test ===")
for criterion_name, score_fn in [
    ("min_covis", lambda t: t['min_covis']),
    ("sum_covis", lambda t: t['sum_covis']),
    ("avg_covis", lambda t: t['avg_covis']),
    ("min*2+avg", lambda t: t['min_covis']*2 + t['avg_covis']),
]:
    # Select triplets
    sorted_tris = sorted(all_triplets, key=score_fn, reverse=True)
    selected_tris = sorted_tris[:111]

    # Extract pairs
    extracted_pairs = extract_pairs_from_triplets(selected_tris)

    # Compare with reference
    matching_pairs = len(extracted_pairs & ref_dense_set)
    tri_matching = sum(1 for t in selected_tris if t['in_ref'])

    print(f"{criterion_name}: triplets={tri_matching}/111, pairs={matching_pairs}/{len(extracted_pairs)} (ref has {len(ref_dense_set)})")

# Try with coverage constraint
print("\n=== Triplet-First with Coverage ===")
for max_img in [2, 3, 4]:
    sorted_tris = sorted(all_triplets, key=lambda t: t['min_covis'], reverse=True)
    selected = []
    img_count = defaultdict(int)

    for t in sorted_tris:
        if len(selected) >= 111:
            break
        tri = t['tri']
        if all(img_count[img] < max_img for img in tri):
            selected.append(t)
            for img in tri:
                img_count[img] += 1

    extracted_pairs = extract_pairs_from_triplets(selected)
    matching_pairs = len(extracted_pairs & ref_dense_set)
    tri_matching = sum(1 for t in selected if t['in_ref'])

    print(f"max_per_img={max_img}: triplets={tri_matching}/{len(selected)}, pairs={matching_pairs}/{len(extracted_pairs)}")
