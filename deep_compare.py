# -*- coding: utf-8 -*-
"""
Deep analysis of why generated pairs differ from original
"""

import os
import re
import math
from collections import defaultdict

def parse_xml_covisibility(xml_path, max_ties=500000):
    """Parse covisibility from XML"""
    covisibility = defaultdict(int)
    tie_count = 0
    current_photo_ids = []
    in_tiepoint = False

    with open(xml_path, 'r', encoding='utf-8') as f:
        for line in f:
            if '<TiePoint>' in line:
                in_tiepoint = True
                current_photo_ids = []
            elif '</TiePoint>' in line:
                in_tiepoint = False
                for i in range(len(current_photo_ids)):
                    for j in range(i+1, len(current_photo_ids)):
                        a, b = current_photo_ids[i], current_photo_ids[j]
                        if a > b:
                            a, b = b, a
                        covisibility[(a, b)] += 1
                tie_count += 1
                if tie_count >= max_ties:
                    break
            elif in_tiepoint and '<PhotoId>' in line:
                m = re.search(r'<PhotoId>(\d+)</PhotoId>', line)
                if m:
                    current_photo_ids.append(int(m.group(1)))
    return covisibility

def parse_pairs_file(path):
    pairs = set()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 1:
                continue
            if len(parts) >= 2:
                try:
                    ids = [int(p) for p in parts]
                    if len(ids) == 2:
                        a, b = min(ids), max(ids)
                        pairs.add((a, b))
                except:
                    pass
    return pairs

def analyze_differences(data_dir, name):
    # Find XML file
    xml_file = None
    for f in os.listdir(data_dir):
        if f.endswith('.xml'):
            xml_file = os.path.join(data_dir, f)
            break

    if not xml_file:
        print(f"XML not found in {data_dir}")
        return

    orig_pairs = parse_pairs_file(os.path.join(data_dir, 'pairs.txt'))
    gen_pairs = parse_pairs_file(os.path.join(data_dir, 'pairs_generated.txt'))

    print(f"\n{'='*60}")
    print(f"Deep Analysis: {name}")
    print(f"{'='*60}")

    # Load covisibility
    print("Loading covisibility data...")
    covis = parse_xml_covisibility(xml_file)
    print(f"Loaded {len(covis)} covisibility pairs")

    # Analyze missed pairs (in original but not in generated)
    missed = orig_pairs - gen_pairs
    extra = gen_pairs - orig_pairs

    # Covisibility analysis for missed vs extra
    missed_covis = [covis.get(p, 0) for p in missed]
    extra_covis = [covis.get(p, 0) for p in extra]
    matching_covis = [covis.get(p, 0) for p in (orig_pairs & gen_pairs)]

    print(f"\nCovisibility Statistics:")
    if missed_covis:
        print(f"  Missed pairs avg covis:   {sum(missed_covis)/len(missed_covis):.1f}")
    if extra_covis:
        print(f"  Extra pairs avg covis:    {sum(extra_covis)/len(extra_covis):.1f}")
    if matching_covis:
        print(f"  Matching pairs avg covis: {sum(matching_covis)/len(matching_covis):.1f}")

    # Build adjacency for triangle analysis
    orig_adj = defaultdict(set)
    gen_adj = defaultdict(set)
    for a, b in orig_pairs:
        orig_adj[a].add(b)
        orig_adj[b].add(a)
    for a, b in gen_pairs:
        gen_adj[a].add(b)
        gen_adj[b].add(a)

    # Check if missed pairs would form triangles in generated graph
    missed_would_form_triangle = 0
    for a, b in missed:
        if gen_adj[a] & gen_adj[b]:  # common neighbors in generated
            missed_would_form_triangle += 1

    # Check if extra pairs form triangles
    extra_forms_triangle = 0
    for a, b in extra:
        if gen_adj[a] & gen_adj[b]:
            extra_forms_triangle += 1

    print(f"\nTriangle Analysis:")
    print(f"  Missed pairs that would form triangles: {missed_would_form_triangle}/{len(missed)} ({100*missed_would_form_triangle/len(missed) if missed else 0:.1f}%)")
    print(f"  Extra pairs that form triangles: {extra_forms_triangle}/{len(extra)} ({100*extra_forms_triangle/len(extra) if extra else 0:.1f}%)")

    # Check 2-hop connectivity of missed pairs in generated graph
    missed_2hop = 0
    for a, b in missed:
        if gen_adj[a] & gen_adj[b]:  # connected via common neighbor
            missed_2hop += 1

    print(f"\n2-Hop Connectivity:")
    print(f"  Missed pairs 2-hop connected in gen: {missed_2hop}/{len(missed)} ({100*missed_2hop/len(missed) if missed else 0:.1f}%)")

    # Degree analysis
    print(f"\nDegree Analysis for nodes with missed pairs:")
    missed_nodes = set()
    for a, b in missed:
        missed_nodes.add(a)
        missed_nodes.add(b)

    gen_degrees = [len(gen_adj[n]) for n in missed_nodes if n in gen_adj]
    orig_degrees = [len(orig_adj[n]) for n in missed_nodes if n in orig_adj]

    if gen_degrees:
        print(f"  Avg degree in generated: {sum(gen_degrees)/len(gen_degrees):.2f}")
    if orig_degrees:
        print(f"  Avg degree in original:  {sum(orig_degrees)/len(orig_degrees):.2f}")

if __name__ == "__main__":
    base_dir = r"D:\codes\cpp\PairSelection\Datas"
    for i in [1]:  # Just analyze first dataset for speed
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            analyze_differences(data_dir, f"Dataset {i}")
