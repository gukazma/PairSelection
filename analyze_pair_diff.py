#!/usr/bin/env python3
"""Analyze differences between reference and generated dense pairs"""

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
print(f"Loaded {len(covis)} pairs")

ref_dense, _, _ = parse_sections('Datas/1/pairs.txt')
gen_dense, _, _ = parse_sections('Datas/1/pairs_generated.txt')

def pair_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

ref_set = set(pair_to_tuple(p) for p in ref_dense)
gen_set = set(pair_to_tuple(p) for p in gen_dense)

matching = ref_set & gen_set
only_ref = ref_set - gen_set
only_gen = gen_set - ref_set

print(f"\nMatching: {len(matching)}, Only ref: {len(only_ref)}, Only gen: {len(only_gen)}")

# Analyze covisibility of different groups
def analyze_group(pairs, name, covis):
    covis_vals = [covis.get(p, 0) for p in pairs]
    id_spans = [p[1] - p[0] for p in pairs]
    print(f"\n{name}:")
    print(f"  Covis - min: {min(covis_vals)}, max: {max(covis_vals)}, avg: {sum(covis_vals)/len(covis_vals):.1f}")
    print(f"  ID span - min: {min(id_spans)}, max: {max(id_spans)}, avg: {sum(id_spans)/len(id_spans):.1f}")

analyze_group(matching, "MATCHING pairs", covis)
analyze_group(only_ref, "Only in REFERENCE", covis)
analyze_group(only_gen, "Only in GENERATED", covis)

# Check if reference pairs are more "local" (smaller ID span)
print("\n=== ID span distribution ===")
ref_spans = sorted([p[1]-p[0] for p in ref_set])
gen_spans = sorted([p[1]-p[0] for p in gen_set])
print(f"Reference - 25th percentile: {ref_spans[len(ref_spans)//4]}, median: {ref_spans[len(ref_spans)//2]}, 75th: {ref_spans[3*len(ref_spans)//4]}")
print(f"Generated - 25th percentile: {gen_spans[len(gen_spans)//4]}, median: {gen_spans[len(gen_spans)//2]}, 75th: {gen_spans[3*len(gen_spans)//4]}")

# Check covisibility distribution
print("\n=== Covisibility distribution ===")
ref_covis = sorted([covis.get(p, 0) for p in ref_set])
gen_covis = sorted([covis.get(p, 0) for p in gen_set])
print(f"Reference - 25th: {ref_covis[len(ref_covis)//4]}, median: {ref_covis[len(ref_covis)//2]}, 75th: {ref_covis[3*len(ref_covis)//4]}")
print(f"Generated - 25th: {gen_covis[len(gen_covis)//4]}, median: {gen_covis[len(gen_covis)//2]}, 75th: {gen_covis[3*len(gen_covis)//4]}")
