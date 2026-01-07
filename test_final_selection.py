#!/usr/bin/env python3
"""
Final optimized selection: Combine per-image rank with intelligent bridge selection.
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

# Build comprehensive pair features
print("Building pair features...")
all_valid_pairs = set()
for p in covis:
    if p[0] not in excluded_images and p[1] not in excluded_images:
        all_valid_pairs.add(p)

pair_features = {}
for p in all_valid_pairs:
    cv = covis[p]
    rank0 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[0]]) if pp == p), 999)
    rank1 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[1]]) if pp == p), 999)
    id_span = abs(p[1] - p[0])

    pair_features[p] = {
        'covis': cv,
        'rank0': rank0,
        'rank1': rank1,
        'min_rank': min(rank0, rank1),
        'max_rank': max(rank0, rank1),
        'sum_rank': rank0 + rank1,
        'id_span': id_span,
        'in_ref': p in ref_dense_set
    }

print(f"Total valid pairs: {len(pair_features)}")

# Current best algorithm (78.8%): per-image rank with multi-pass
def current_best(pair_features, max_degree=4, target=198):
    # Sort by min_rank, then covis
    sorted_pairs = sorted(pair_features.items(),
                          key=lambda x: (x[1]['min_rank'], -x[1]['covis']))

    img_degree = defaultdict(int)
    selected = set()

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

    return selected

selected = current_best(pair_features)
matching = len(selected & ref_dense_set)
print(f"\nCurrent best: {matching}/198 = {matching*100/198:.1f}%")

# Analyze what we're missing
missed = ref_dense_set - selected
extra = selected - ref_dense_set

print(f"Missed: {len(missed)}, Extra: {len(extra)}")

# Key insight: For images with degree constraint hit, we might be selecting
# the wrong pairs. Let's see if we can swap some extra pairs for missed pairs.
print("\n=== Swap Analysis ===")

# For each missed pair, can we swap it with an extra pair?
def can_swap(missed_pair, extra_pair, selected, img_degree, max_degree=4):
    """Check if swapping missed for extra is valid"""
    # After removing extra
    new_deg = img_degree.copy()
    new_deg[extra_pair[0]] -= 1
    new_deg[extra_pair[1]] -= 1

    # Can we add missed?
    if new_deg[missed_pair[0]] >= max_degree or new_deg[missed_pair[1]] >= max_degree:
        return False
    return True

# Rebuild selection with degree tracking
img_degree = defaultdict(int)
selected = set()
sorted_pairs = sorted(pair_features.items(), key=lambda x: (x[1]['min_rank'], -x[1]['covis']))

for p, f in sorted_pairs:
    if len(selected) >= 198:
        break
    if img_degree[p[0]] == 0 and img_degree[p[1]] == 0:
        selected.add(p)
        img_degree[p[0]] += 1
        img_degree[p[1]] += 1

for p, f in sorted_pairs:
    if len(selected) >= 198:
        break
    if p in selected:
        continue
    if img_degree[p[0]] == 0 or img_degree[p[1]] == 0:
        if img_degree[p[0]] < 4 and img_degree[p[1]] < 4:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

for p, f in sorted_pairs:
    if len(selected) >= 198:
        break
    if p in selected:
        continue
    if img_degree[p[0]] < 4 and img_degree[p[1]] < 4:
        selected.add(p)
        img_degree[p[0]] += 1
        img_degree[p[1]] += 1

missed = ref_dense_set - selected
extra = selected - ref_dense_set

# Find possible swaps
possible_swaps = []
for m in missed:
    m_f = pair_features[m]
    for e in extra:
        e_f = pair_features[e]
        # Would swapping improve? Check if missed is "better" than extra
        # Better = lower min_rank or same min_rank but higher covis
        if can_swap(m, e, selected, img_degree):
            # Is the missed pair better connected to the images?
            score_improvement = (e_f['min_rank'] - m_f['min_rank']) * 100 + (m_f['covis'] - e_f['covis']) / 10
            if score_improvement > 0:
                possible_swaps.append((m, e, score_improvement))

possible_swaps.sort(key=lambda x: -x[2])
print(f"Possible beneficial swaps: {len(possible_swaps)}")
if possible_swaps:
    print("Top 10 swaps:")
    for m, e, score in possible_swaps[:10]:
        m_f = pair_features[m]
        e_f = pair_features[e]
        print(f"  Swap {e} (ranks={e_f['min_rank']},{e_f['max_rank']}, cv={e_f['covis']})")
        print(f"  For  {m} (ranks={m_f['min_rank']},{m_f['max_rank']}, cv={m_f['covis']})")

# Try greedy swapping
print("\n=== Greedy Swap ===")
for m, e, score in possible_swaps:
    if e in selected and m not in selected:
        # Perform swap
        selected.remove(e)
        selected.add(m)
        img_degree[e[0]] -= 1
        img_degree[e[1]] -= 1
        img_degree[m[0]] += 1
        img_degree[m[1]] += 1

matching = len(selected & ref_dense_set)
print(f"After swaps: {matching}/198 = {matching*100/198:.1f}%")

# Alternative approach: Iterative refinement
print("\n=== Iterative Refinement ===")

def iterative_selection(pair_features, max_degree=4, target=198, iterations=10):
    # Start with current best
    sorted_pairs = sorted(pair_features.items(),
                          key=lambda x: (x[1]['min_rank'], -x[1]['covis']))

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

    matching = len(selected & ref_dense_set)
    print(f"  Initial: {matching}/198")

    # Iterative improvement
    for iter in range(iterations):
        improved = False
        missed = ref_dense_set - selected
        extra = selected - ref_dense_set

        for m in missed:
            for e in extra:
                # Check if swap is valid
                new_deg = img_degree.copy()
                new_deg[e[0]] -= 1
                new_deg[e[1]] -= 1
                if new_deg[m[0]] < max_degree and new_deg[m[1]] < max_degree:
                    # Do the swap
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

    matching = len(selected & ref_dense_set)
    print(f"  After {iter+1} iterations: {matching}/198")
    return selected

selected = iterative_selection(pair_features)
matching = len(selected & ref_dense_set)
print(f"\nIterative: {matching}/198 = {matching*100/198:.1f}%")

# Final analysis
print("\n=== Final Analysis ===")
missed = ref_dense_set - selected
extra = selected - ref_dense_set
print(f"Still missing: {len(missed)}")

# Why are these pairs missed?
print("\nMissed pairs characteristics:")
for p in list(missed)[:20]:
    f = pair_features[p]
    # What images are involved?
    deg0 = img_degree[p[0]]
    deg1 = img_degree[p[1]]
    print(f"  {p}: ranks=({f['rank0']},{f['rank1']}), covis={f['covis']}, img_degrees=({deg0},{deg1})")
