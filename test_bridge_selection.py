#!/usr/bin/env python3
"""
Bridge-aware pair selection.
Key insight: Reference includes bridge pairs (large id_span) even with lower covis.
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

# Analyze bridge pairs in reference
print("=== Bridge Pairs in Reference ===")
bridge_threshold = 500  # ID span threshold for "bridge" pairs
ref_bridges = []
ref_local = []
for p in ref_dense_set:
    id_span = abs(p[1] - p[0])
    if id_span > bridge_threshold:
        ref_bridges.append((p, id_span, covis[p]))
    else:
        ref_local.append((p, id_span, covis[p]))

print(f"Bridge pairs (span > {bridge_threshold}): {len(ref_bridges)}")
print(f"Local pairs: {len(ref_local)}")

ref_bridges.sort(key=lambda x: -x[1])
print(f"\nTop bridge pairs:")
for p, span, cv in ref_bridges[:15]:
    rank0 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[0]]) if pp == p), 999)
    rank1 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[1]]) if pp == p), 999)
    print(f"  {p}: span={span}, covis={cv}, ranks=({rank0},{rank1})")

# Build pair features
K = 6
candidate_pairs = set()
for img in valid_images:
    pairs = sorted_pairs_map[img]
    for p, cv in pairs[:K]:
        candidate_pairs.add(p)

pair_features = {}
for p in candidate_pairs:
    cv = covis[p]
    rank0 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[0]]) if pp == p), 999)
    rank1 = next((i for i, (pp, _) in enumerate(sorted_pairs_map[p[1]]) if pp == p), 999)
    id_span = abs(p[1] - p[0])
    is_bridge = id_span > bridge_threshold

    pair_features[p] = {
        'covis': cv,
        'rank0': rank0,
        'rank1': rank1,
        'min_rank': min(rank0, rank1),
        'max_rank': max(rank0, rank1),
        'id_span': id_span,
        'is_bridge': is_bridge,
        'in_ref': p in ref_dense_set
    }

# How many bridge pairs are in candidates?
candidate_bridges = [p for p in pair_features if pair_features[p]['is_bridge']]
ref_bridges_in_cand = [p for p in candidate_bridges if pair_features[p]['in_ref']]
print(f"\nBridge pairs in candidates: {len(candidate_bridges)}")
print(f"Reference bridges in candidates: {len(ref_bridges_in_cand)}")

# Bridge-aware selection
print("\n=== Bridge-Aware Selection ===")

def bridge_aware_selection(pair_features, max_degree=4, target=198, bridge_priority=True):
    """
    Phase 1: Select bridge pairs first
    Phase 2: Fill with local pairs
    """
    img_degree = defaultdict(int)
    selected = set()

    bridges = [(p, f) for p, f in pair_features.items() if f['is_bridge']]
    local = [(p, f) for p, f in pair_features.items() if not f['is_bridge']]

    # Sort bridges by min_rank, then covis
    bridges.sort(key=lambda x: (x[1]['min_rank'], -x[1]['covis']))
    local.sort(key=lambda x: (x[1]['min_rank'], -x[1]['covis']))

    if bridge_priority:
        # Phase 1: Select all valid bridges
        for p, f in bridges:
            if len(selected) >= target:
                break
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

        print(f"  After bridge phase: {len(selected)} pairs")

    # Phase 2: Multi-pass local selection
    for p, f in local:
        if len(selected) >= target:
            break
        if img_degree[p[0]] == 0 and img_degree[p[1]] == 0:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    for p, f in local:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] == 0 or img_degree[p[1]] == 0:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

    for p, f in local:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    # If still not enough, add any remaining
    all_sorted = sorted(pair_features.items(), key=lambda x: (x[1]['min_rank'], -x[1]['covis']))
    for p, f in all_sorted:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    return selected

# Test with bridge priority
selected = bridge_aware_selection(pair_features, bridge_priority=True)
matching = len(selected & ref_dense_set)
print(f"Bridge-aware: {matching}/198 = {matching*100/198:.1f}%")

# Test without bridge priority (baseline)
selected = bridge_aware_selection(pair_features, bridge_priority=False)
matching = len(selected & ref_dense_set)
print(f"Without bridge priority: {matching}/198 = {matching*100/198:.1f}%")

# Interleaved bridge selection
print("\n=== Interleaved Bridge Selection ===")

def interleaved_selection(pair_features, max_degree=4, target=198):
    """
    Interleave bridge and local pairs based on coverage needs.
    """
    img_degree = defaultdict(int)
    selected = set()

    # Separate bridges and local
    bridges = [(p, f) for p, f in pair_features.items() if f['is_bridge']]
    local = [(p, f) for p, f in pair_features.items() if not f['is_bridge']]

    bridges.sort(key=lambda x: (x[1]['min_rank'], -x[1]['covis']))
    local.sort(key=lambda x: (x[1]['min_rank'], -x[1]['covis']))

    # Phase 1: Both degree 0 - mix of bridges and local
    all_pairs = bridges + local
    all_pairs.sort(key=lambda x: (x[1]['min_rank'], -x[1]['covis']))

    for p, f in all_pairs:
        if len(selected) >= target:
            break
        if img_degree[p[0]] == 0 and img_degree[p[1]] == 0:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    print(f"  After phase 1 (both deg 0): {len(selected)} pairs")

    # Phase 2: Prioritize bridges where one image is uncovered
    uncovered = valid_images - set(img_degree.keys())
    for p, f in bridges:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if p[0] in uncovered or p[1] in uncovered:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

    print(f"  After phase 2 (bridge coverage): {len(selected)} pairs")

    # Phase 3: One degree 0
    for p, f in all_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] == 0 or img_degree[p[1]] == 0:
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1

    # Phase 4: Fill
    for p, f in all_pairs:
        if len(selected) >= target:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    return selected

selected = interleaved_selection(pair_features)
matching = len(selected & ref_dense_set)
print(f"Interleaved: {matching}/198 = {matching*100/198:.1f}%")

# Analyze which reference bridges we're missing
print("\n=== Analyzing Missed Reference Bridges ===")
selected_bridges = [p for p in selected if pair_features[p]['is_bridge']]
ref_bridges_set = set(p for p, _, _ in ref_bridges)
matching_bridges = set(selected_bridges) & ref_bridges_set
missed_bridges = ref_bridges_set - set(selected_bridges)
extra_bridges = set(selected_bridges) - ref_bridges_set

print(f"Selected bridges: {len(selected_bridges)}")
print(f"Reference bridges: {len(ref_bridges_set)}")
print(f"Matching bridges: {len(matching_bridges)}")
print(f"Missed bridges: {len(missed_bridges)}")
print(f"Extra bridges: {len(extra_bridges)}")

if missed_bridges:
    print("\nMissed reference bridges:")
    for p in list(missed_bridges)[:10]:
        f = pair_features[p]
        print(f"  {p}: span={f['id_span']}, covis={f['covis']}, ranks=({f['rank0']},{f['rank1']})")

if extra_bridges:
    print("\nExtra non-ref bridges:")
    for p in list(extra_bridges)[:10]:
        f = pair_features[p]
        print(f"  {p}: span={f['id_span']}, covis={f['covis']}, ranks=({f['rank0']},{f['rank1']})")

# Try varying bridge threshold
print("\n=== Testing Bridge Thresholds ===")
for threshold in [300, 400, 500, 600, 800, 1000]:
    # Rebuild features with new threshold
    for p in pair_features:
        pair_features[p]['is_bridge'] = pair_features[p]['id_span'] > threshold

    selected = bridge_aware_selection(pair_features)
    matching = len(selected & ref_dense_set)
    n_bridges = sum(1 for p in selected if pair_features[p]['is_bridge'])
    print(f"  threshold={threshold}: {matching}/198 = {matching*100/198:.1f}%, bridges={n_bridges}")
