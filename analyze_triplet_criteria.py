#!/usr/bin/env python3
"""
Deep analysis of reference triplet selection criteria.
What makes a triplet get selected in the reference?
"""

import re
from collections import defaultdict
import glob

def parse_xml_covisibility(xml_path):
    """Parse XML to get covisibility counts"""
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

print("Loading XML...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])
print(f"Loaded {len(covis)} covisibility pairs")

ref_dense, _, ref_triplets = parse_sections('Datas/1/pairs.txt')
gen_dense, _, gen_triplets = parse_sections('Datas/1/pairs_generated.txt')

def triplet_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

def pair_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

ref_tri_set = set(triplet_to_tuple(t) for t in ref_triplets)
gen_tri_set = set(triplet_to_tuple(t) for t in gen_triplets)
ref_dense_set = set(pair_to_tuple(p) for p in ref_dense)
gen_dense_set = set(pair_to_tuple(p) for p in gen_dense)

# For each image, count how many times it appears in triplets
ref_img_count = defaultdict(int)
for t in ref_tri_set:
    for img in t:
        ref_img_count[img] += 1

gen_img_count = defaultdict(int)
for t in gen_tri_set:
    for img in t:
        gen_img_count[img] += 1

print("\n=== Image coverage in triplets ===")
print(f"Reference: {len(ref_img_count)} images, min count: {min(ref_img_count.values())}, max: {max(ref_img_count.values())}, avg: {sum(ref_img_count.values())/len(ref_img_count):.1f}")
print(f"Generated: {len(gen_img_count)} images, min count: {min(gen_img_count.values())}, max: {max(gen_img_count.values())}, avg: {sum(gen_img_count.values())/len(gen_img_count):.1f}")

# Images that appear only once in reference triplets - these are "edge" images
ref_once = [img for img, cnt in ref_img_count.items() if cnt == 1]
gen_once = [img for img, cnt in gen_img_count.items() if cnt == 1]
print(f"\nImages appearing once in triplets - Ref: {len(ref_once)}, Gen: {len(gen_once)}")

# Analyze triplets by their "missing edge" covisibility
# The edge NOT in dense pairs - what's its covisibility?
def get_missing_edge_covis(t, dense_set, covis):
    a, b, c = t
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    for e in edges:
        if e not in dense_set:
            return covis.get(e, 0)
    return -1  # All edges in dense (shouldn't happen)

print("\n=== Missing edge analysis ===")
ref_missing = [get_missing_edge_covis(t, ref_dense_set, covis) for t in ref_tri_set]
print(f"Reference missing edge covis - min: {min(ref_missing)}, max: {max(ref_missing)}, avg: {sum(ref_missing)/len(ref_missing):.1f}")

# For generated, we might have different dense pairs structure
gen_missing = [get_missing_edge_covis(t, gen_dense_set, covis) for t in gen_tri_set]
print(f"Generated missing edge covis - min: {min(gen_missing)}, max: {max(gen_missing)}, avg: {sum(gen_missing)/len(gen_missing):.1f}")

# Analyze the "kept" edges in reference
def get_kept_edges_covis(t, dense_set, covis):
    a, b, c = t
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    kept = [covis.get(e, 0) for e in edges if e in dense_set]
    return kept

print("\n=== Kept edges analysis ===")
all_kept_ref = []
for t in ref_tri_set:
    all_kept_ref.extend(get_kept_edges_covis(t, ref_dense_set, covis))
print(f"Reference kept edge covis - min: {min(all_kept_ref)}, max: {max(all_kept_ref)}, avg: {sum(all_kept_ref)/len(all_kept_ref):.1f}")

# Is the missing edge always the one with lowest covisibility?
def is_missing_lowest(t, dense_set, covis):
    a, b, c = t
    edges = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]
    covis_values = [covis.get(e, 0) for e in edges]
    min_covis = min(covis_values)
    min_idx = covis_values.index(min_covis)
    missing_edge = [e for e in edges if e not in dense_set]
    if not missing_edge:
        return None
    return edges[min_idx] == missing_edge[0]

print("\n=== Is missing edge always the lowest covis? ===")
results = [is_missing_lowest(t, ref_dense_set, covis) for t in ref_tri_set]
true_count = sum(1 for r in results if r == True)
false_count = sum(1 for r in results if r == False)
none_count = sum(1 for r in results if r is None)
print(f"Yes: {true_count}, No: {false_count}, N/A: {none_count}")
print(f"Missing edge is lowest covis: {true_count/(true_count+false_count)*100:.1f}%")
