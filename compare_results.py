# -*- coding: utf-8 -*-
"""
Compare generated pairs with original reference pairs
"""

import os
from collections import defaultdict

def parse_pairs_file(path):
    """Parse pairs.txt format and return set of pairs"""
    pairs = set()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 1:
                continue  # photo count line
            if len(parts) >= 2:
                try:
                    ids = [int(p) for p in parts]
                    if len(ids) == 2:
                        a, b = min(ids), max(ids)
                        pairs.add((a, b))
                    else:
                        for i in range(len(ids)):
                            for j in range(i+1, len(ids)):
                                a, b = min(ids[i], ids[j]), max(ids[i], ids[j])
                                pairs.add((a, b))
                except:
                    pass
    return pairs

def analyze_dataset(data_dir, name):
    orig_path = os.path.join(data_dir, 'pairs.txt')
    gen_path = os.path.join(data_dir, 'pairs_generated.txt')

    if not os.path.exists(orig_path) or not os.path.exists(gen_path):
        print(f"Files not found in {data_dir}")
        return

    orig_pairs = parse_pairs_file(orig_path)
    gen_pairs = parse_pairs_file(gen_path)

    matching = orig_pairs & gen_pairs
    only_orig = orig_pairs - gen_pairs
    only_gen = gen_pairs - orig_pairs

    precision = len(matching) / len(gen_pairs) * 100 if gen_pairs else 0
    recall = len(matching) / len(orig_pairs) * 100 if orig_pairs else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n{'='*60}")
    print(f"Dataset: {name}")
    print(f"{'='*60}")
    print(f"Original pairs:   {len(orig_pairs)}")
    print(f"Generated pairs:  {len(gen_pairs)}")
    print(f"Matching pairs:   {len(matching)}")
    print(f"Only in original: {len(only_orig)} (missed)")
    print(f"Only in generated: {len(only_gen)} (extra)")
    print(f"\nPrecision: {precision:.2f}%")
    print(f"Recall:    {recall:.2f}%")
    print(f"F1 Score:  {f1:.2f}%")

    # Degree distribution comparison
    orig_degree = defaultdict(int)
    gen_degree = defaultdict(int)
    for a, b in orig_pairs:
        orig_degree[a] += 1
        orig_degree[b] += 1
    for a, b in gen_pairs:
        gen_degree[a] += 1
        gen_degree[b] += 1

    orig_avg_degree = sum(orig_degree.values()) / len(orig_degree) if orig_degree else 0
    gen_avg_degree = sum(gen_degree.values()) / len(gen_degree) if gen_degree else 0

    print(f"\nOriginal avg degree: {orig_avg_degree:.2f}")
    print(f"Generated avg degree: {gen_avg_degree:.2f}")

    # Show sample differences
    print(f"\nSample pairs only in original (first 5):")
    for p in sorted(list(only_orig))[:5]:
        print(f"  {p[0]} - {p[1]}")

    print(f"\nSample pairs only in generated (first 5):")
    for p in sorted(list(only_gen))[:5]:
        print(f"  {p[0]} - {p[1]}")

if __name__ == "__main__":
    base_dir = r"D:\codes\cpp\PairSelection\Datas"
    for i in [1, 2, 3, 4]:
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            analyze_dataset(data_dir, f"Dataset {i}")
