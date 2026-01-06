# -*- coding: utf-8 -*-
"""
Comprehensive comparison analysis
"""

import os
import re
from collections import defaultdict

def parse_xml_photos(xml_path, max_photos=5000):
    photos = {}
    current_photo = {}
    current_id = None
    in_photo = False
    in_pose = False
    in_rotation = False
    in_center = False
    photo_count = 0

    with open(xml_path, 'r', encoding='utf-8') as f:
        for line in f:
            if '<Photo>' in line:
                in_photo = True
                current_photo = {'rotation': {}, 'center': {}}
            elif '</Photo>' in line:
                if current_id is not None and 'x' in current_photo['center']:
                    photos[current_id] = current_photo
                    photo_count += 1
                    if photo_count >= max_photos:
                        break
                in_photo = False
                current_id = None
            elif in_photo:
                if '<Id>' in line:
                    m = re.search(r'<Id>(\d+)</Id>', line)
                    if m:
                        current_id = int(m.group(1))
                elif '<Pose>' in line:
                    in_pose = True
                elif '</Pose>' in line:
                    in_pose = False
                elif in_pose:
                    if '<Rotation>' in line:
                        in_rotation = True
                    elif '</Rotation>' in line:
                        in_rotation = False
                    elif '<Center>' in line:
                        in_center = True
                    elif '</Center>' in line:
                        in_center = False
                    elif in_rotation:
                        for key in ['M_02', 'M_12', 'M_22']:
                            if f'<{key}>' in line:
                                m = re.search(rf'<{key}>([^<]+)</{key}>', line)
                                if m:
                                    current_photo['rotation'][key] = float(m.group(1))
                    elif in_center:
                        for key in ['x', 'y', 'z']:
                            if f'<{key}>' in line:
                                m = re.search(rf'<{key}>([^<]+)</{key}>', line)
                                if m:
                                    current_photo['center'][key] = float(m.group(1))
                if '<MedianDepth>' in line:
                    m = re.search(r'<MedianDepth>([^<]+)</MedianDepth>', line)
                    if m:
                        current_photo['median_depth'] = float(m.group(1))
    return photos

def parse_xml_covisibility(xml_path, max_ties=500000):
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

def calc_angle(v1, v2):
    import math
    dot = sum(a*b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a*a for a in v1))
    mag2 = math.sqrt(sum(b*b for b in v2))
    if mag1 == 0 or mag2 == 0:
        return 0
    cos_angle = max(-1, min(1, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))

def get_optical_axis(photo):
    R = photo['rotation']
    if 'M_02' not in R:
        return None
    return [R['M_02'], R['M_12'], R['M_22']]

def comprehensive_analysis(data_dir, name):
    xml_file = None
    for f in os.listdir(data_dir):
        if f.endswith('.xml'):
            xml_file = os.path.join(data_dir, f)
            break

    if not xml_file:
        return

    print(f"\n{'='*70}")
    print(f"Comprehensive Analysis: {name}")
    print(f"{'='*70}")

    # Load data
    photos = parse_xml_photos(xml_file)
    covis = parse_xml_covisibility(xml_file)
    orig_pairs = parse_pairs_file(os.path.join(data_dir, 'pairs.txt'))
    gen_pairs = parse_pairs_file(os.path.join(data_dir, 'pairs_generated.txt'))

    print(f"Photos: {len(photos)}, Covis pairs: {len(covis)}")
    print(f"Original: {len(orig_pairs)}, Generated: {len(gen_pairs)}")

    # Filter to valid pairs (both photos exist in XML)
    valid_photo_ids = set(photos.keys())
    orig_valid = {p for p in orig_pairs if p[0] in valid_photo_ids and p[1] in valid_photo_ids}
    gen_valid = {p for p in gen_pairs if p[0] in valid_photo_ids and p[1] in valid_photo_ids}

    print(f"Valid original: {len(orig_valid)}, Valid generated: {len(gen_valid)}")

    matching = orig_valid & gen_valid
    missed = orig_valid - gen_valid
    extra = gen_valid - orig_valid

    print(f"Matching: {len(matching)}, Missed: {len(missed)}, Extra: {len(extra)}")

    # Analyze characteristics
    def analyze_pair_set(pair_set, label):
        if not pair_set:
            return
        covis_vals = [covis.get(p, 0) for p in pair_set]
        angles = []
        for a, b in pair_set:
            if a in photos and b in photos:
                axis1 = get_optical_axis(photos[a])
                axis2 = get_optical_axis(photos[b])
                if axis1 and axis2:
                    angles.append(calc_angle(axis1, axis2))

        print(f"\n{label}:")
        print(f"  Count: {len(pair_set)}")
        print(f"  Avg covisibility: {sum(covis_vals)/len(covis_vals):.1f}")
        print(f"  Max covisibility: {max(covis_vals)}")
        print(f"  Min covisibility: {min(covis_vals)}")
        if angles:
            print(f"  Avg convergence angle: {sum(angles)/len(angles):.2f}°")
            print(f"  Max convergence angle: {max(angles):.2f}°")

    analyze_pair_set(matching, "Matching pairs (both agree)")
    analyze_pair_set(missed, "Missed pairs (in orig, not in gen)")
    analyze_pair_set(extra, "Extra pairs (in gen, not in orig)")

    # Rank analysis - where do missed pairs rank?
    print(f"\nRank Analysis for Missed Pairs:")
    all_neighbors = defaultdict(list)
    for pair, count in covis.items():
        a, b = pair
        all_neighbors[a].append((b, count))
        all_neighbors[b].append((a, count))

    for node in all_neighbors:
        all_neighbors[node].sort(key=lambda x: -x[1])

    missed_ranks = []
    for a, b in missed:
        # Find rank of b in a's neighbor list
        for rank, (neighbor, _) in enumerate(all_neighbors[a]):
            if neighbor == b:
                missed_ranks.append(rank + 1)
                break

    if missed_ranks:
        print(f"  Avg rank of missed pairs: {sum(missed_ranks)/len(missed_ranks):.1f}")
        print(f"  Max rank: {max(missed_ranks)}")
        print(f"  Pairs in top-5: {sum(1 for r in missed_ranks if r <= 5)} ({100*sum(1 for r in missed_ranks if r <= 5)/len(missed_ranks):.1f}%)")

if __name__ == "__main__":
    base_dir = r"D:\codes\cpp\PairSelection\Datas"
    for i in [1, 2, 3, 4]:
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            comprehensive_analysis(data_dir, f"Dataset {i}")
