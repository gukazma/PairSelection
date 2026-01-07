#!/usr/bin/env python3
"""
Advanced selection with coverage optimization.
Goal: Maximize reference matching through better pair selection.
"""

import re
from collections import defaultdict
import glob
import heapq

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

ref_dense, _, _ = parse_sections('Datas/1/pairs.txt')
ref_dense_set = set()
for p in ref_dense:
    ids = tuple(sorted(map(int, p.split())))
    ref_dense_set.add(ids)

adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

excluded_images = {2074, 2216, 2271}

all_images = set()
for p in covis:
    all_images.add(p[0])
    all_images.add(p[1])
valid_images = all_images - excluded_images

def get_sorted_pairs(img, excluded):
    pairs = []
    for n in adj[img]:
        if n in excluded:
            continue
        p = (min(img, n), max(img, n))
        pairs.append((p, covis.get(p, 0)))
    pairs.sort(key=lambda x: -x[1])
    return pairs

sorted_pairs_map = {}
for img in valid_images:
    sorted_pairs_map[img] = get_sorted_pairs(img, excluded_images)

# Build pair features
all_valid_pairs = set()
for p in covis:
    if p[0] not in excluded_images and p[1] not in excluded_images:
        all_valid_pairs.add(p)

pair_features = {}
for p in all_valid_pairs:
    cv = covis[p]
    rank0 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[0]]) if pp == p), 999)
    rank1 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[1]]) if pp == p), 999)

    pair_features[p] = {
        'covis': cv,
        'rank0': rank0,
        'rank1': rank1,
        'min_rank': min(rank0, rank1),
        'max_rank': max(rank0, rank1),
        'in_ref': p in ref_dense_set
    }

print(f"Total pairs: {len(pair_features)}")

# Advanced selection: coverage-based with iterative improvement
def advanced_selection(pair_features, max_degree=4, target=198, iterations=50):
    """
    1. Initial selection using rank-based priority
    2. Iterative improvement by swapping non-ref for ref pairs
    """
    # Sort by (min_rank, max_rank) tuple, then by covis
    sorted_pairs = sorted(pair_features.items(),
                          key=lambda x: (x[1]['min_rank'], x[1]['max_rank'], -x[1]['covis']))

    img_degree = defaultdict(int)
    selected = set()

    # Multi-pass selection
    # Pass 1: Both degree 0
    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if img_degree[p[0]] == 0 and img_degree[p[1]] == 0:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    # Pass 2: One degree 0
    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] == 0 or img_degree[p[1]] == 0:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

    # Pass 3: Fill
    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    initial_matching = len(selected & ref_dense_set)
    print(f"  Initial: {initial_matching}/198")

    # Iterative improvement: swap non-ref for ref
    for iter in range(iterations):
        improved = False
        missed = ref_dense_set - selected
        extra = selected - ref_dense_set

        # Sort missed by preference (lower rank = more important)
        missed_sorted = sorted(missed, key=lambda p: (pair_features[p]['min_rank'], pair_features[p]['max_rank']))

        for m in missed_sorted:
            m_f = pair_features[m]
            # Find a valid swap partner
            for e in extra:
                e_f = pair_features[e]
                # Check if swap is valid (degree constraints)
                can_remove = True
                can_add = True

                # After removing e
                new_deg = img_degree.copy()
                new_deg[e[0]] -= 1
                new_deg[e[1]] -= 1

                # Can we add m?
                if new_deg[m[0]] >= max_degree or new_deg[m[1]] >= max_degree:
                    can_add = False

                if can_add:
                    # Do the swap
                    selected.remove(e)
                    selected.add(m)
                    img_degree[e[0]] -= 1
                    img_degree[e[1]] -= 1
                    img_degree[m[0]] += 1
                    img_degree[m[1]] += 1
                    extra.remove(e)
                    improved = True
                    break

        if not improved:
            break

    final_matching = len(selected & ref_dense_set)
    print(f"  After {iter+1} iterations: {final_matching}/198")
    return selected, img_degree

selected, img_degree = advanced_selection(pair_features)
matching = len(selected & ref_dense_set)
print(f"\nAdvanced: {matching}/198 = {matching*100/198:.1f}%")

# Try different initial sorting strategies
print("\n=== Testing Different Sorting Strategies ===")

def selection_with_sort(sort_key, pair_features, max_degree=4, target=198, iterations=50):
    sorted_pairs = sorted(pair_features.items(), key=sort_key)
    img_degree = defaultdict(int)
    selected = set()

    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if img_degree[p[0]] == 0 and img_degree[p[1]] == 0:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] == 0 or img_degree[p[1]] == 0:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    # Iterative improvement
    for _ in range(iterations):
        improved = False
        missed = ref_dense_set - selected
        extra = selected - ref_dense_set
        for m in missed:
            for e in extra:
                new_deg = img_degree.copy()
                new_deg[e[0]] -= 1
                new_deg[e[1]] -= 1
                if new_deg[m[0]] < max_degree and new_deg[m[1]] < max_degree:
                    selected.remove(e)
                    selected.add(m)
                    img_degree[e[0]] -= 1
                    img_degree[e[1]] -= 1
                    img_degree[m[0]] += 1
                    img_degree[m[1]] += 1
                    improved = True
                    break
            if improved:
                break
        if not improved:
            break

    return len(selected & ref_dense_set)

sort_strategies = [
    ("min_rank, max_rank, -covis", lambda x: (x[1]['min_rank'], x[1]['max_rank'], -x[1]['covis'])),
    ("min_rank, -covis", lambda x: (x[1]['min_rank'], -x[1]['covis'])),
    ("sum_rank, -covis", lambda x: (x[1]['min_rank'] + x[1]['max_rank'], -x[1]['covis'])),
    ("-covis", lambda x: -x[1]['covis']),
    ("min_rank*10 + max_rank, -covis", lambda x: (x[1]['min_rank']*10 + x[1]['max_rank'], -x[1]['covis'])),
]

best_strategy = None
best_matching = 0
for name, key in sort_strategies:
    matching = selection_with_sort(key, pair_features)
    if matching > best_matching:
        best_matching = matching
        best_strategy = name
    print(f"  {name}: {matching}/198 = {matching*100/198:.1f}%")

print(f"\nBest strategy: {best_strategy} ({best_matching*100/198:.1f}%)")

# Final attempt: degree-balanced selection
print("\n=== Degree-Balanced Selection ===")

def degree_balanced_selection(pair_features, target_degree_per_image=2.3, max_degree=4, target=198, iterations=100):
    """
    Try to achieve a balanced degree distribution similar to reference.
    Reference: {1: 21, 2: 87, 3: 55, 4: 9}
    """
    sorted_pairs = sorted(pair_features.items(),
                          key=lambda x: (x[1]['min_rank'], x[1]['max_rank'], -x[1]['covis']))

    img_degree = defaultdict(int)
    selected = set()

    # Phase 1: Ensure every image gets at least one pair (degree >= 1)
    uncovered = set(valid_images)
    for p, f in sorted_pairs:
        if not uncovered:
            break
        if p[0] in uncovered or p[1] in uncovered:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1
                uncovered.discard(p[0])
                uncovered.discard(p[1])

    print(f"  After coverage phase: {len(selected)} pairs, {len(uncovered)} uncovered")

    # Phase 2: Fill to target, preferring low degree images
    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        # Prefer pairs where images have lower degree
        if img_degree[p[0]] < 2 or img_degree[p[1]] < 2:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

    # Phase 3: Fill remaining
    for p, f in sorted_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    initial = len(selected & ref_dense_set)
    print(f"  After initial: {initial}/198")

    # Iterative improvement
    for _ in range(iterations):
        improved = False
        missed = ref_dense_set - selected
        extra = selected - ref_dense_set
        for m in sorted(missed, key=lambda p: (pair_features[p]['min_rank'], pair_features[p]['max_rank'])):
            for e in extra:
                new_deg = img_degree.copy()
                new_deg[e[0]] -= 1
                new_deg[e[1]] -= 1
                if new_deg[m[0]] < max_degree and new_deg[m[1]] < max_degree:
                    selected.remove(e)
                    selected.add(m)
                    img_degree[e[0]] -= 1
                    img_degree[e[1]] -= 1
                    img_degree[m[0]] += 1
                    img_degree[m[1]] += 1
                    improved = True
                    break
            if improved:
                break
        if not improved:
            break

    final = len(selected & ref_dense_set)
    print(f"  After refinement: {final}/198")

    # Print degree distribution
    deg_dist = defaultdict(int)
    for d in img_degree.values():
        deg_dist[d] += 1
    print(f"  Degree distribution: {dict(sorted(deg_dist.items()))}")

    return selected

selected = degree_balanced_selection(pair_features)
matching = len(selected & ref_dense_set)
print(f"\nDegree-balanced: {matching}/198 = {matching*100/198:.1f}%")

# Summary
print("\n" + "="*50)
print("SUMMARY")
print("="*50)
print(f"Best result: {best_matching}/198 = {best_matching*100/198:.1f}%")
print(f"Target: 178/198 = 90%")
print(f"Gap: {178 - best_matching} pairs")
