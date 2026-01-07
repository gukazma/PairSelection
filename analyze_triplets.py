#!/usr/bin/env python3
"""Deep analysis of reference triplets vs generated triplets"""
import xml.etree.ElementTree as ET
import re
from collections import defaultdict
import math

def parse_xml_covisibility(xml_path):
    """Parse XML to get covisibility counts"""
    covis = defaultdict(int)

    with open(xml_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Find all TiePoints
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
    """Parse pairs.txt format"""
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
            if section_idx == 0:
                dense = section_data
            elif section_idx == 1:
                refine = section_data
            elif section_idx == 2:
                triplets = section_data
            section_idx += 1
            i += count + 1
        else:
            i += 1
    return dense, refine, triplets

print("Loading XML covisibility (this may take a moment)...")
import glob
xml_files = glob.glob('Datas/1/*.xml')
covis = parse_xml_covisibility(xml_files[0]) if xml_files else {}
print(f"Loaded {len(covis)} covisibility pairs")

# Parse reference
ref_dense, ref_refine, ref_triplets = parse_sections('Datas/1/pairs.txt')
gen_dense, gen_refine, gen_triplets = parse_sections('Datas/1/pairs_generated.txt')

def triplet_to_tuple(line):
    return tuple(sorted(map(int, line.split())))

def get_triplet_edges(t):
    """Get sorted edges of a triplet"""
    a, b, c = t
    return [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]

def analyze_triplet(t, covis):
    """Get statistics for a triplet"""
    edges = get_triplet_edges(t)
    covis_counts = [covis.get(e, 0) for e in edges]
    return {
        'min_covis': min(covis_counts),
        'max_covis': max(covis_counts),
        'avg_covis': sum(covis_counts) / 3,
        'covis_counts': covis_counts,
        'id_span': max(t) - min(t)
    }

ref_tri_set = set(triplet_to_tuple(t) for t in ref_triplets)
gen_tri_set = set(triplet_to_tuple(t) for t in gen_triplets)

matching = ref_tri_set & gen_tri_set
only_ref = ref_tri_set - gen_tri_set
only_gen = gen_tri_set - ref_tri_set

print(f"\n=== Reference triplet analysis ===")
ref_stats = [analyze_triplet(t, covis) for t in ref_tri_set]
print(f"Min covis - min: {min(s['min_covis'] for s in ref_stats)}, max: {max(s['min_covis'] for s in ref_stats)}, avg: {sum(s['min_covis'] for s in ref_stats)/len(ref_stats):.1f}")
print(f"Avg covis - min: {min(s['avg_covis'] for s in ref_stats):.1f}, max: {max(s['avg_covis'] for s in ref_stats):.1f}, avg: {sum(s['avg_covis'] for s in ref_stats)/len(ref_stats):.1f}")
print(f"ID span - min: {min(s['id_span'] for s in ref_stats)}, max: {max(s['id_span'] for s in ref_stats)}, avg: {sum(s['id_span'] for s in ref_stats)/len(ref_stats):.1f}")

print(f"\n=== Generated triplet analysis ===")
gen_stats = [analyze_triplet(t, covis) for t in gen_tri_set]
print(f"Min covis - min: {min(s['min_covis'] for s in gen_stats)}, max: {max(s['min_covis'] for s in gen_stats)}, avg: {sum(s['min_covis'] for s in gen_stats)/len(gen_stats):.1f}")
print(f"Avg covis - min: {min(s['avg_covis'] for s in gen_stats):.1f}, max: {max(s['avg_covis'] for s in gen_stats):.1f}, avg: {sum(s['avg_covis'] for s in gen_stats)/len(gen_stats):.1f}")
print(f"ID span - min: {min(s['id_span'] for s in gen_stats)}, max: {max(s['id_span'] for s in gen_stats)}, avg: {sum(s['id_span'] for s in gen_stats)/len(gen_stats):.1f}")

print(f"\n=== Triplets ONLY in reference ===")
only_ref_stats = [analyze_triplet(t, covis) for t in only_ref]
print(f"Min covis - avg: {sum(s['min_covis'] for s in only_ref_stats)/len(only_ref_stats):.1f}")
print(f"Avg covis - avg: {sum(s['avg_covis'] for s in only_ref_stats)/len(only_ref_stats):.1f}")
print(f"ID span - avg: {sum(s['id_span'] for s in only_ref_stats)/len(only_ref_stats):.1f}")

print(f"\n=== Triplets ONLY in generated ===")
only_gen_stats = [analyze_triplet(t, covis) for t in only_gen]
print(f"Min covis - avg: {sum(s['min_covis'] for s in only_gen_stats)/len(only_gen_stats):.1f}")
print(f"Avg covis - avg: {sum(s['avg_covis'] for s in only_gen_stats)/len(only_gen_stats):.1f}")
print(f"ID span - avg: {sum(s['id_span'] for s in only_gen_stats)/len(only_gen_stats):.1f}")

# Check image coverage
print("\n=== Image coverage analysis ===")
ref_images = set()
for t in ref_tri_set:
    ref_images.update(t)
gen_images = set()
for t in gen_tri_set:
    gen_images.update(t)
print(f"Reference covers {len(ref_images)} unique images")
print(f"Generated covers {len(gen_images)} unique images")
print(f"Image overlap: {len(ref_images & gen_images)}")
