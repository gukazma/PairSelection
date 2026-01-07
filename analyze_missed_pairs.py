#!/usr/bin/env python3
"""
Deep analysis of missed pairs to find the missing criterion.
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
ref_dense_set = set()
for p in ref_dense:
    ids = tuple(sorted(map(int, p.split())))
    ref_dense_set.add(ids)

ref_triplet_set = set()
for t in ref_triplets:
    ids = tuple(sorted(map(int, t.split())))
    ref_triplet_set.add(ids)

# Build adjacency
adj = defaultdict(set)
for p in covis:
    adj[p[0]].add(p[1])
    adj[p[1]].add(p[0])

excluded_images = {2074, 2216, 2271}

# Get rank for a pair relative to an image
def get_rank(img, pair, excluded):
    neighbors = [(n, covis[(min(img,n), max(img,n))]) for n in adj[img] if n not in excluded]
    neighbors.sort(key=lambda x: -x[1])
    neighbor_pairs = [(min(img, n), max(img, n)) for n, _ in neighbors]
    try:
        return neighbor_pairs.index(pair)
    except ValueError:
        return 999

# Analyze reference pairs by their rank profile
print("=== Reference Pair Rank Analysis ===")
ref_rank_profiles = []
for p in ref_dense_set:
    rank1 = get_rank(p[0], p, excluded_images)
    rank2 = get_rank(p[1], p, excluded_images)
    min_rank = min(rank1, rank2)
    max_rank = max(rank1, rank2)
    cv = covis.get(p, 0)
    ref_rank_profiles.append((p, rank1, rank2, min_rank, max_rank, cv))

# Sort by max rank
ref_rank_profiles.sort(key=lambda x: x[4])
print("\nReference pairs by max rank:")
rank_groups = defaultdict(list)
for p, r1, r2, min_r, max_r, cv in ref_rank_profiles:
    rank_groups[max_r].append((p, r1, r2, cv))

for max_r in sorted(rank_groups.keys())[:10]:
    count = len(rank_groups[max_r])
    sample = rank_groups[max_r][:3]
    print(f"  max_rank={max_r}: {count} pairs, e.g., {sample}")

# Key insight: What's the maximum max_rank in reference?
max_max_rank = max(x[4] for x in ref_rank_profiles)
print(f"\nMaximum max_rank in reference: {max_max_rank}")

# How many pairs have max_rank > 5?
high_rank = [x for x in ref_rank_profiles if x[4] > 5]
print(f"Pairs with max_rank > 5: {len(high_rank)}")
for p, r1, r2, min_r, max_r, cv in high_rank[:10]:
    print(f"  {p}: ranks=({r1},{r2}), covis={cv}")

# Hypothesis: Reference requires min_rank <= threshold
print("\n=== Testing min_rank constraint ===")
for threshold in [0, 1, 2, 3, 4, 5]:
    qualifying = [x for x in ref_rank_profiles if x[3] <= threshold]
    print(f"  min_rank <= {threshold}: {len(qualifying)}/198 pairs qualify")

# Analyze triplet connection
print("\n=== Pairs in Triplets ===")
# For each reference pair, check if it's part of a triplet
pair_in_triplet = defaultdict(set)
for tri in ref_triplet_set:
    a, b, c = tri
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    for e in edges:
        if e in ref_dense_set:
            pair_in_triplet[e].add(tri)

pairs_with_triplet = [p for p in ref_dense_set if p in pair_in_triplet]
print(f"Reference pairs that are in at least one triplet: {len(pairs_with_triplet)}/198")

# What about pairs NOT in any triplet?
pairs_no_triplet = [p for p in ref_dense_set if p not in pair_in_triplet]
print(f"Reference pairs NOT in any triplet: {len(pairs_no_triplet)}")
for p in pairs_no_triplet[:10]:
    rank1 = get_rank(p[0], p, excluded_images)
    rank2 = get_rank(p[1], p, excluded_images)
    print(f"  {p}: ranks=({rank1},{rank2}), covis={covis.get(p,0)}")

# Analyze image degree in reference
print("\n=== Image Degree Analysis ===")
ref_adj = defaultdict(set)
for p in ref_dense_set:
    ref_adj[p[0]].add(p[1])
    ref_adj[p[1]].add(p[0])

degrees = [(img, len(ref_adj[img])) for img in ref_adj]
degrees.sort(key=lambda x: x[1])

degree_dist = defaultdict(int)
for img, d in degrees:
    degree_dist[d] += 1
print(f"Degree distribution: {dict(sorted(degree_dist.items()))}")

# Check: for each image, are selected pairs always top-k by covis?
print("\n=== Per-image selection pattern ===")
selection_ranks = []
for img, neighbors in ref_adj.items():
    # Get all neighbors sorted by covis
    all_neighbors = sorted(
        [(n, covis[(min(img,n), max(img,n))]) for n in adj[img] if n not in excluded_images],
        key=lambda x: -x[1]
    )
    # Find ranks of selected neighbors
    ranks = []
    for n in neighbors:
        for idx, (nn, cv) in enumerate(all_neighbors):
            if nn == n:
                ranks.append(idx)
                break
    selection_ranks.extend(ranks)

    if len(ranks) <= 3:
        print(f"  Img {img}: degree={len(neighbors)}, ranks={sorted(ranks)}")

print(f"\nAll selection ranks: min={min(selection_ranks)}, max={max(selection_ranks)}, avg={sum(selection_ranks)/len(selection_ranks):.1f}")

# Distribution
sel_rank_dist = defaultdict(int)
for r in selection_ranks:
    sel_rank_dist[r] += 1
print(f"Selection rank distribution: {dict(sorted(sel_rank_dist.items())[:15])}")

# Key hypothesis: For each image, select its top-N neighbors (N = target degree)
print("\n=== Per-Image Top-N Selection ===")
def select_per_image_topn(excluded, degree_per_image=2, max_pairs=198):
    """For each image, greedily select its top neighbors up to degree limit"""
    img_degree = defaultdict(int)
    selected = set()

    # Collect all candidate pairs with their "importance"
    candidates = []
    for img in adj:
        if img in excluded:
            continue
        # Get image's top neighbors
        neighbors = sorted(
            [(n, covis[(min(img,n), max(img,n))]) for n in adj[img] if n not in excluded],
            key=lambda x: -x[1]
        )
        for rank, (n, cv) in enumerate(neighbors):
            p = (min(img, n), max(img, n))
            # Score: prioritize being a top choice for both images
            candidates.append((p, rank, img, n, cv))

    # Sort by rank (lower is better), then by covis
    candidates.sort(key=lambda x: (x[1], -x[4]))

    # Remove duplicates but keep track of which ranks they appear at
    pair_ranks = defaultdict(list)
    for p, rank, img, n, cv in candidates:
        pair_ranks[p].append((rank, cv))

    # Score pairs by min rank
    pair_scores = []
    for p, ranks in pair_ranks.items():
        min_rank = min(r for r, cv in ranks)
        max_rank = max(r for r, cv in ranks)
        cv = ranks[0][1]
        pair_scores.append((p, min_rank, max_rank, cv))

    # Sort by min_rank, then max_rank, then covis
    pair_scores.sort(key=lambda x: (x[1], x[2], -x[3]))

    # Select with degree constraint
    for p, min_r, max_r, cv in pair_scores:
        if len(selected) >= max_pairs:
            break
        if img_degree[p[0]] < degree_per_image and img_degree[p[1]] < degree_per_image:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    # Phase 2: relax degree constraint
    for p, min_r, max_r, cv in pair_scores:
        if len(selected) >= max_pairs:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < 4 and img_degree[p[1]] < 4:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    return selected

for deg in [2, 3, 4]:
    selected = select_per_image_topn(excluded_images, degree_per_image=deg)
    matching = len(selected & ref_dense_set)
    print(f"  degree_per_image={deg}: selected={len(selected)}, matching={matching}/198={matching*100/198:.1f}%")

# Alternative: ensure each image has at least one pair in its top-2
print("\n=== Guaranteed Top-2 Coverage ===")
def select_with_coverage_guarantee(excluded, max_pairs=198, max_degree=4):
    """
    Phase 1: Ensure each image has at least one pair from its top-2
    Phase 2: Fill remaining with best available
    """
    img_degree = defaultdict(int)
    selected = set()
    covered = set()  # Images that have at least one pair from top-2

    all_images = set()
    for p in covis:
        if p[0] not in excluded:
            all_images.add(p[0])
        if p[1] not in excluded:
            all_images.add(p[1])

    # Get top-2 pairs for each image
    top2_map = {}
    for img in all_images:
        neighbors = sorted(
            [(n, covis[(min(img,n), max(img,n))]) for n in adj[img] if n not in excluded],
            key=lambda x: -x[1]
        )[:2]
        top2_map[img] = [(min(img, n), max(img, n)) for n, cv in neighbors]

    # Phase 1: Cover all images with their top-2 pairs
    # Prioritize pairs that are top-2 for both images
    all_pairs = []
    for p in covis:
        if p[0] in excluded or p[1] in excluded:
            continue
        is_top2_for_0 = p in top2_map.get(p[0], [])
        is_top2_for_1 = p in top2_map.get(p[1], [])
        all_pairs.append((p, covis[p], is_top2_for_0, is_top2_for_1))

    # Sort: both top-2 first, then one top-2, then by covis
    all_pairs.sort(key=lambda x: (-(x[2] and x[3]), -(x[2] or x[3]), -x[1]))

    for p, cv, t0, t1 in all_pairs:
        if len(selected) >= max_pairs:
            break
        if p in selected:
            continue

        # Check if this pair helps cover uncovered images
        helps_0 = p[0] not in covered and t0
        helps_1 = p[1] not in covered and t1

        if helps_0 or helps_1 or (t0 and t1):
            if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
                selected.add(p)
                img_degree[p[0]] += 1
                img_degree[p[1]] += 1
                if t0:
                    covered.add(p[0])
                if t1:
                    covered.add(p[1])

    # Phase 2: Fill remaining
    all_pairs.sort(key=lambda x: -x[1])  # Sort by covis
    for p, cv, t0, t1 in all_pairs:
        if len(selected) >= max_pairs:
            break
        if p in selected:
            continue
        if img_degree[p[0]] < max_degree and img_degree[p[1]] < max_degree:
            selected.add(p)
            img_degree[p[0]] += 1
            img_degree[p[1]] += 1

    return selected, covered

selected, covered = select_with_coverage_guarantee(excluded_images)
matching = len(selected & ref_dense_set)
print(f"With coverage guarantee: selected={len(selected)}, covered={len(covered)} images, matching={matching}/198={matching*100/198:.1f}%")
