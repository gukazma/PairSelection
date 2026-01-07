#!/usr/bin/env python3
"""
Analyze what makes a pair get selected in reference dense pairs.
Build a model to predict reference pair selection.
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

print("Loading data...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])

ref_dense, _, _ = parse_sections('Datas/1/pairs.txt')
gen_dense, _, _ = parse_sections('Datas/1/pairs_generated.txt')

ref_set = set(pair_to_tuple(p) for p in ref_dense)
gen_set = set(pair_to_tuple(p) for p in gen_dense)

# Build adjacency
adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

# Collect all images
all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])

# For each pair, calculate various features
print("Calculating pair features...")
pair_features = []
for pair, cv in covis.items():
    a, b = pair
    id_span = b - a

    # How many common neighbors?
    common_neighbors = len(adj[a] & adj[b])

    # Degree of each image
    deg_a = len(adj[a])
    deg_b = len(adj[b])

    # Is this pair part of a triangle with high covis?
    max_tri_covis = 0
    min_tri_covis = float('inf')
    for c in adj[a] & adj[b]:
        e1 = (min(a,c), max(a,c))
        e2 = (min(b,c), max(b,c))
        tri_min = min(cv, covis.get(e1, 0), covis.get(e2, 0))
        tri_max = max(cv, covis.get(e1, 0), covis.get(e2, 0))
        max_tri_covis = max(max_tri_covis, tri_min)  # Best weakest link
        if tri_min > 0:
            min_tri_covis = min(min_tri_covis, tri_min)

    if min_tri_covis == float('inf'):
        min_tri_covis = 0

    in_ref = pair in ref_set

    pair_features.append({
        'pair': pair,
        'covis': cv,
        'id_span': id_span,
        'common_neighbors': common_neighbors,
        'deg_a': deg_a,
        'deg_b': deg_b,
        'max_tri_covis': max_tri_covis,
        'min_tri_covis': min_tri_covis,
        'in_ref': in_ref
    })

# Analyze correlations
print("\n=== Feature analysis ===")
ref_pairs = [p for p in pair_features if p['in_ref']]
non_ref_pairs = [p for p in pair_features if not p['in_ref']]

print(f"Reference pairs: {len(ref_pairs)}")
print(f"Non-reference pairs: {len(non_ref_pairs)}")

def stats(pairs, feature):
    vals = [p[feature] for p in pairs]
    return min(vals), max(vals), sum(vals)/len(vals)

for feature in ['covis', 'id_span', 'common_neighbors', 'deg_a', 'max_tri_covis']:
    ref_stats = stats(ref_pairs, feature)
    non_stats = stats(non_ref_pairs, feature)
    print(f"\n{feature}:")
    print(f"  Reference: min={ref_stats[0]}, max={ref_stats[1]}, avg={ref_stats[2]:.1f}")
    print(f"  Non-ref:   min={non_stats[0]}, max={non_stats[1]}, avg={non_stats[2]:.1f}")

# Try to find a selection rule that matches reference
print("\n=== Testing selection rules ===")

# Excluded images
excluded_images = {2074, 2216, 2271}

# Rule 1: Exclude problematic images, then select by covis
valid_pairs = [p for p in pair_features if p['pair'][0] not in excluded_images and p['pair'][1] not in excluded_images]
valid_pairs.sort(key=lambda x: -x['covis'])
selected = set(p['pair'] for p in valid_pairs[:198])
matching = len(selected & ref_set)
print(f"Rule 1 (exclude 3 imgs, top 198 by covis): {matching}/198 = {matching*100/198:.1f}%")

# Rule 2: Same but with degree constraint
img_degree = defaultdict(int)
selected = set()
for p in valid_pairs:
    pair = p['pair']
    if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
        selected.add(pair)
        img_degree[pair[0]] += 1
        img_degree[pair[1]] += 1
        if len(selected) >= 198:
            break
matching = len(selected & ref_set)
print(f"Rule 2 (+ degree constraint): {matching}/198 = {matching*100/198:.1f}%")

# Rule 3: Prioritize pairs that are part of good triangles
valid_pairs.sort(key=lambda x: -x['max_tri_covis'])
img_degree = defaultdict(int)
selected = set()
for p in valid_pairs:
    pair = p['pair']
    if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
        selected.add(pair)
        img_degree[pair[0]] += 1
        img_degree[pair[1]] += 1
        if len(selected) >= 198:
            break
matching = len(selected & ref_set)
print(f"Rule 3 (by best triangle quality): {matching}/198 = {matching*100/198:.1f}%")

# Rule 4: Combined score
for w_covis, w_tri in [(0.5, 0.5), (0.3, 0.7), (0.7, 0.3)]:
    valid_pairs.sort(key=lambda x: -(w_covis * x['covis'] + w_tri * x['max_tri_covis']))
    img_degree = defaultdict(int)
    selected = set()
    for p in valid_pairs:
        pair = p['pair']
        if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
            selected.add(pair)
            img_degree[pair[0]] += 1
            img_degree[pair[1]] += 1
            if len(selected) >= 198:
                break
    matching = len(selected & ref_set)
    print(f"Rule 4 (w_covis={w_covis}, w_tri={w_tri}): {matching}/198 = {matching*100/198:.1f}%")

# Rule 5: Multi-pass selection (prioritize coverage)
print("\n=== Multi-pass with excluded images ===")
img_degree = defaultdict(int)
selected = set()
valid_pairs_by_covis = sorted([p for p in valid_pairs], key=lambda x: -x['covis'])

# Pass 1: Pairs where both images have degree 0
for p in valid_pairs_by_covis:
    if len(selected) >= 198: break
    pair = p['pair']
    if img_degree[pair[0]] == 0 and img_degree[pair[1]] == 0:
        if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
            selected.add(pair)
            img_degree[pair[0]] += 1
            img_degree[pair[1]] += 1

# Pass 2: Pairs where at least one image has degree 0
for p in valid_pairs_by_covis:
    if len(selected) >= 198: break
    pair = p['pair']
    if pair in selected: continue
    if img_degree[pair[0]] == 0 or img_degree[pair[1]] == 0:
        if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
            selected.add(pair)
            img_degree[pair[0]] += 1
            img_degree[pair[1]] += 1

# Pass 3: Fill remaining
for p in valid_pairs_by_covis:
    if len(selected) >= 198: break
    pair = p['pair']
    if pair in selected: continue
    if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
        selected.add(pair)
        img_degree[pair[0]] += 1
        img_degree[pair[1]] += 1

matching = len(selected & ref_set)
print(f"Multi-pass (exclude 3 imgs): {matching}/198 = {matching*100/198:.1f}%")

# Check which pairs are in reference but we missed
missed = ref_set - selected
print(f"\nMissed {len(missed)} reference pairs")
missed_features = [p for p in pair_features if p['pair'] in missed]
missed_features.sort(key=lambda x: -x['covis'])
print("Top 10 missed (by covis):")
for p in missed_features[:10]:
    print(f"  {p['pair']}: covis={p['covis']}, id_span={p['id_span']}, common_neighbors={p['common_neighbors']}")

# What pairs did we select that aren't in reference?
extra = selected - ref_set
print(f"\nSelected {len(extra)} pairs NOT in reference")
extra_features = [p for p in pair_features if p['pair'] in extra]
extra_features.sort(key=lambda x: -x['covis'])
print("Top 10 extra (by covis):")
for p in extra_features[:10]:
    print(f"  {p['pair']}: covis={p['covis']}, id_span={p['id_span']}, common_neighbors={p['common_neighbors']}")
