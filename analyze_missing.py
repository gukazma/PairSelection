#!/usr/bin/env python3
"""Analyze specific missing reference pairs"""

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

print("Loading XML...")
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0])

ref_dense, _, ref_triplets = parse_sections('Datas/1/pairs.txt')
gen_dense, _, _ = parse_sections('Datas/1/pairs_generated.txt')

def pair_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

ref_set = set(pair_to_tuple(p) for p in ref_dense)
gen_set = set(pair_to_tuple(p) for p in gen_dense)

only_ref = sorted(ref_set - gen_set)

print(f"\nPairs ONLY in reference ({len(only_ref)}):")
print("Pair\t\tCovis\tID_Span")
for p in only_ref[:20]:
    c = covis.get(p, 0)
    span = p[1] - p[0]
    print(f"{p[0]}-{p[1]}\t\t{c}\t{span}")

# Check if these pairs are in candidate list (have covisibility)
print(f"\nChecking if missing pairs are valid candidates...")
missing_not_candidate = []
for p in only_ref:
    if covis.get(p, 0) == 0:
        missing_not_candidate.append(p)
print(f"Missing pairs with zero covisibility: {len(missing_not_candidate)}")

# Check score ranking of missing pairs
all_candidates = [(covis.get(p, 0), p) for p in ref_set | gen_set if covis.get(p, 0) > 0]
all_candidates.sort(reverse=True)

print(f"\nRanking of missing pairs among all candidates ({len(all_candidates)} total):")
for p in only_ref[:10]:
    c = covis.get(p, 0)
    rank = next((i for i, (cv, pr) in enumerate(all_candidates) if pr == p), -1)
    print(f"{p}: covis={c}, rank={rank+1}/{len(all_candidates)}")
