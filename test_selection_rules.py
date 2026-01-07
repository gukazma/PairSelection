#!/usr/bin/env python3
"""
Test selection rules that prioritize common neighbors.
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
ref_set = set(pair_to_tuple(p) for p in ref_dense)

# Build adjacency
adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

# Excluded images
excluded_images = {2074, 2216, 2271}

# Calculate features for all pairs
pair_data = []
max_covis = max(covis.values())
for pair, cv in covis.items():
    if pair[0] in excluded_images or pair[1] in excluded_images:
        continue
    common = len(adj[pair[0]] & adj[pair[1]])
    id_span = pair[1] - pair[0]
    pair_data.append({
        'pair': pair,
        'covis': cv,
        'common': common,
        'id_span': id_span
    })

print(f"Valid pairs: {len(pair_data)}")

def test_rule(pair_data, score_fn, name, ref_set, max_degree=4):
    # Sort by score
    sorted_pairs = sorted(pair_data, key=score_fn, reverse=True)

    img_degree = defaultdict(int)
    selected = set()

    # Multi-pass selection
    for threshold in [0, 0, 1, 999]:
        for p in sorted_pairs:
            if len(selected) >= 198: break
            pair = p['pair']
            if pair in selected: continue

            if threshold == 999:
                # Final pass: any valid pair
                if img_degree[pair[0]] < max_degree and img_degree[pair[1]] < max_degree:
                    selected.add(pair)
                    img_degree[pair[0]] += 1
                    img_degree[pair[1]] += 1
            else:
                if threshold == 0 and img_degree[pair[0]] == 0 and img_degree[pair[1]] == 0:
                    if img_degree[pair[0]] < max_degree and img_degree[pair[1]] < max_degree:
                        selected.add(pair)
                        img_degree[pair[0]] += 1
                        img_degree[pair[1]] += 1
                elif threshold == 1 and (img_degree[pair[0]] == 0 or img_degree[pair[1]] == 0):
                    if img_degree[pair[0]] < max_degree and img_degree[pair[1]] < max_degree:
                        selected.add(pair)
                        img_degree[pair[0]] += 1
                        img_degree[pair[1]] += 1

    matching = len(selected & ref_set)
    precision = matching / len(selected) if selected else 0
    recall = matching / len(ref_set) if ref_set else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    print(f"{name}: {matching}/198 = {matching*100/198:.1f}%, P={precision*100:.1f}%, R={recall*100:.1f}%, F1={f1*100:.1f}%")
    return selected

# Test various scoring functions
print("\n=== Testing score functions ===")

# Pure covisibility
test_rule(pair_data, lambda p: p['covis'], "covis only", ref_set)

# Common neighbors only
test_rule(pair_data, lambda p: p['common'], "common_neighbors only", ref_set)

# Combined scores
test_rule(pair_data, lambda p: p['covis'] * p['common'], "covis * common", ref_set)
test_rule(pair_data, lambda p: p['covis'] + p['common'] * 10, "covis + 10*common", ref_set)
test_rule(pair_data, lambda p: p['covis'] + p['common'] * 20, "covis + 20*common", ref_set)
test_rule(pair_data, lambda p: p['covis'] + p['common'] * 50, "covis + 50*common", ref_set)

# Normalized
max_common = max(p['common'] for p in pair_data)
test_rule(pair_data, lambda p: 0.5 * p['covis']/max_covis + 0.5 * p['common']/max_common, "0.5*norm_covis + 0.5*norm_common", ref_set)
test_rule(pair_data, lambda p: 0.3 * p['covis']/max_covis + 0.7 * p['common']/max_common, "0.3*norm_covis + 0.7*norm_common", ref_set)
test_rule(pair_data, lambda p: 0.7 * p['covis']/max_covis + 0.3 * p['common']/max_common, "0.7*norm_covis + 0.3*norm_common", ref_set)

# With ID span penalty
test_rule(pair_data, lambda p: p['covis'] - p['id_span'] * 0.1, "covis - 0.1*id_span", ref_set)
test_rule(pair_data, lambda p: p['covis'] * p['common'] / (1 + p['id_span']/100), "covis*common / (1+span/100)", ref_set)

# Triangle-based: max covis of any triangle containing this edge
print("\n=== Triangle-aware scoring ===")
# Precompute best triangle covis for each pair
for p in pair_data:
    pair = p['pair']
    a, b = pair
    best_tri = 0
    for c in adj[a] & adj[b]:
        e1 = (min(a,c), max(a,c))
        e2 = (min(b,c), max(b,c))
        tri_min = min(p['covis'], covis.get(e1, 0), covis.get(e2, 0))
        best_tri = max(best_tri, tri_min)
    p['best_tri'] = best_tri

test_rule(pair_data, lambda p: p['best_tri'], "best_triangle_min_edge", ref_set)
test_rule(pair_data, lambda p: p['covis'] + p['best_tri'], "covis + best_tri", ref_set)
test_rule(pair_data, lambda p: p['common'] * p['best_tri'], "common * best_tri", ref_set)

# Find the best combination
print("\n=== Grid search for optimal weights ===")
best_f1 = 0
best_params = None
for w_covis in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:
    for w_common in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        for w_tri in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:
            if w_covis + w_common + w_tri == 0:
                continue
            def score_fn(p, wc=w_covis, wcm=w_common, wt=w_tri):
                return wc * p['covis']/max_covis + wcm * p['common']/max_common + wt * p['best_tri']/max_covis

            sorted_pairs = sorted(pair_data, key=score_fn, reverse=True)
            img_degree = defaultdict(int)
            selected = set()

            for threshold in [0, 1, 999]:
                for pa in sorted_pairs:
                    if len(selected) >= 198: break
                    pair = pa['pair']
                    if pair in selected: continue
                    if threshold == 999:
                        if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
                            selected.add(pair)
                            img_degree[pair[0]] += 1
                            img_degree[pair[1]] += 1
                    elif threshold == 0 and img_degree[pair[0]] == 0 and img_degree[pair[1]] == 0:
                        if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
                            selected.add(pair)
                            img_degree[pair[0]] += 1
                            img_degree[pair[1]] += 1
                    elif threshold == 1 and (img_degree[pair[0]] == 0 or img_degree[pair[1]] == 0):
                        if img_degree[pair[0]] < 4 and img_degree[pair[1]] < 4:
                            selected.add(pair)
                            img_degree[pair[0]] += 1
                            img_degree[pair[1]] += 1

            matching = len(selected & ref_set)
            precision = matching / len(selected) if selected else 0
            recall = matching / len(ref_set)
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

            if f1 > best_f1:
                best_f1 = f1
                best_params = (w_covis, w_common, w_tri, matching, precision, recall)

print(f"Best: w_covis={best_params[0]}, w_common={best_params[1]}, w_tri={best_params[2]}")
print(f"  Matching: {best_params[3]}/198, P={best_params[4]*100:.1f}%, R={best_params[5]*100:.1f}%, F1={best_f1*100:.1f}%")
