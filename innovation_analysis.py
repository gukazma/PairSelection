# -*- coding: utf-8 -*-
"""
深入分析选择逻辑的创新点
1. 非严格Top-K的原因
2. 冗余避免机制
3. 连通性保证机制
"""

import re
import math
import os
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
    photogroup_id = 0

    with open(xml_path, 'r', encoding='utf-8') as f:
        for line in f:
            if '<Photogroup>' in line or '<PhotoGroup>' in line:
                photogroup_id += 1
            if '<Photo>' in line:
                in_photo = True
                current_photo = {'rotation': {}, 'center': {}, 'photogroup': photogroup_id}
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
                        for key in ['M_00', 'M_01', 'M_02', 'M_10', 'M_11', 'M_12', 'M_20', 'M_21', 'M_22']:
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

def get_optical_axis(photo):
    R = photo['rotation']
    if 'M_20' not in R:
        return None
    return [R['M_02'], R['M_12'], R['M_22']]

def calc_angle(v1, v2):
    dot = sum(a*b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a*a for a in v1))
    mag2 = math.sqrt(sum(b*b for b in v2))
    if mag1 == 0 or mag2 == 0:
        return 0
    cos_angle = max(-1, min(1, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))

def calc_convergence_angle(photo1, photo2):
    axis1 = get_optical_axis(photo1)
    axis2 = get_optical_axis(photo2)
    if axis1 is None or axis2 is None:
        return None
    return calc_angle(axis1, axis2)

def calc_gsd_ratio(photo1, photo2):
    if 'median_depth' not in photo1 or 'median_depth' not in photo2:
        return None
    d1, d2 = photo1['median_depth'], photo2['median_depth']
    if d1 == 0 or d2 == 0:
        return None
    return max(d1, d2) / min(d1, d2)

def parse_pairs(pairs_path):
    pairs = []
    with open(pairs_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    ids = [int(p) for p in parts]
                    for i in range(len(ids)):
                        for j in range(i+1, len(ids)):
                            pairs.append((ids[i], ids[j]))
                except:
                    pass
    return pairs

def parse_covisibility(xml_path, max_ties=100000):
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

def analyze_innovation(data_dir, name):
    """分析选择逻辑的创新点"""
    print(f"\n{'='*70}")
    print(f"创新点深度分析: {name}")
    print('='*70)

    xml_file = None
    pairs_file = None
    for f in os.listdir(data_dir):
        if f.endswith('.xml'):
            xml_file = os.path.join(data_dir, f)
        if f == 'pairs.txt':
            pairs_file = os.path.join(data_dir, f)

    if not xml_file or not pairs_file:
        return

    photos = parse_xml_photos(xml_file)
    selected_pairs = parse_pairs(pairs_file)
    selected_set = set((min(a,b), max(a,b)) for a, b in selected_pairs)
    covisibility = parse_covisibility(xml_file)

    print(f"照片数: {len(photos)}, 选中像对: {len(selected_set)}, 共视像对: {len(covisibility)}")

    # 构建图
    all_neighbors = defaultdict(list)
    for pair, covis in covisibility.items():
        a, b = pair
        if a not in photos or b not in photos:
            continue
        p1, p2 = photos[a], photos[b]
        angle = calc_convergence_angle(p1, p2)
        gsd = calc_gsd_ratio(p1, p2)
        is_selected = pair in selected_set

        all_neighbors[a].append({'id': b, 'covis': covis, 'angle': angle, 'gsd': gsd, 'selected': is_selected})
        all_neighbors[b].append({'id': a, 'covis': covis, 'angle': angle, 'gsd': gsd, 'selected': is_selected})

    # 创新点1: 分析"非Top-K"选择的原因
    print("\n" + "="*50)
    print("创新点1: 智能非Top-K选择分析")
    print("="*50)

    skip_topk_reasons = defaultdict(int)

    for photo_id, neighbors in all_neighbors.items():
        if len(neighbors) < 5:
            continue

        sorted_neighbors = sorted(neighbors, key=lambda x: -x['covis'])
        selected_ids = set(n['id'] for n in neighbors if n['selected'])

        for rank, n in enumerate(sorted_neighbors[:10]):  # 检查前10名
            if n['id'] not in selected_ids:  # 高排名但未选中
                # 分析原因
                if n['angle'] is not None and n['angle'] > 15:
                    skip_topk_reasons['交会角过大(>15°)'] += 1
                elif n['gsd'] is not None and n['gsd'] > 1.3:
                    skip_topk_reasons['GSD比过大(>1.3)'] += 1
                else:
                    # 检查是否是冗余连接
                    # 检查n['id']和photo_id是否已经通过其他路径连接
                    common_selected = selected_ids & set(nn['id'] for nn in all_neighbors[n['id']] if nn['selected'])
                    if common_selected:
                        skip_topk_reasons['存在共同邻居(冗余)'] += 1
                    else:
                        skip_topk_reasons['其他原因'] += 1

    print("高排名但未被选中的原因分布:")
    total = sum(skip_topk_reasons.values())
    for reason, count in sorted(skip_topk_reasons.items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count} ({100*count/total:.1f}%)")

    # 创新点2: 三角形闭合分析
    print("\n" + "="*50)
    print("创新点2: 三角形闭合模式分析")
    print("="*50)

    # 统计选中图中的三角形数量
    selected_graph = defaultdict(set)
    for pair in selected_set:
        a, b = pair
        selected_graph[a].add(b)
        selected_graph[b].add(a)

    triangle_count = 0
    triangle_nodes = set()
    for a in selected_graph:
        for b in selected_graph[a]:
            if b > a:
                common = selected_graph[a] & selected_graph[b]
                for c in common:
                    if c > b:
                        triangle_count += 1
                        triangle_nodes.update([a, b, c])

    print(f"  选中图中的三角形数量: {triangle_count}")
    print(f"  参与三角形的节点数: {len(triangle_nodes)} ({100*len(triangle_nodes)/len(selected_graph):.1f}%)")
    print(f"  平均每个三角形对应的节点数: {3*triangle_count/max(len(triangle_nodes),1):.2f}")

    # 创新点3: 链式vs星型结构分析
    print("\n" + "="*50)
    print("创新点3: 局部拓扑结构分析")
    print("="*50)

    # 计算聚类系数
    clustering_coeffs = []
    for node in selected_graph:
        neighbors = list(selected_graph[node])
        if len(neighbors) < 2:
            continue
        # 计算邻居之间的连接数
        edges_between_neighbors = 0
        for i in range(len(neighbors)):
            for j in range(i+1, len(neighbors)):
                if neighbors[j] in selected_graph[neighbors[i]]:
                    edges_between_neighbors += 1
        possible_edges = len(neighbors) * (len(neighbors) - 1) / 2
        cc = edges_between_neighbors / possible_edges if possible_edges > 0 else 0
        clustering_coeffs.append(cc)

    if clustering_coeffs:
        avg_cc = sum(clustering_coeffs) / len(clustering_coeffs)
        print(f"  平均聚类系数: {avg_cc:.3f}")
        print(f"  聚类系数范围: [{min(clustering_coeffs):.3f}, {max(clustering_coeffs):.3f}]")

        # 聚类系数分布
        cc_dist = {'<0.1': 0, '0.1-0.3': 0, '0.3-0.5': 0, '>0.5': 0}
        for cc in clustering_coeffs:
            if cc < 0.1:
                cc_dist['<0.1'] += 1
            elif cc < 0.3:
                cc_dist['0.1-0.3'] += 1
            elif cc < 0.5:
                cc_dist['0.3-0.5'] += 1
            else:
                cc_dist['>0.5'] += 1
        print(f"  聚类系数分布:")
        for k, v in cc_dist.items():
            print(f"    {k}: {v} ({100*v/len(clustering_coeffs):.1f}%)")

    # 创新点4: 综合评分函数反推
    print("\n" + "="*50)
    print("创新点4: 评分函数反推")
    print("="*50)

    # 对于相同重叠度区间的像对，分析选中vs未选中的差异
    covis_ranges = [(50, 100), (100, 200), (200, 300)]

    for low, high in covis_ranges:
        selected_in_range = []
        not_selected_in_range = []

        for pair, covis in covisibility.items():
            if not (low <= covis < high):
                continue
            a, b = pair
            if a not in photos or b not in photos:
                continue

            p1, p2 = photos[a], photos[b]
            angle = calc_convergence_angle(p1, p2)
            gsd = calc_gsd_ratio(p1, p2)

            if angle is None or gsd is None:
                continue

            features = {'covis': covis, 'angle': angle, 'gsd': gsd}

            if pair in selected_set:
                selected_in_range.append(features)
            else:
                not_selected_in_range.append(features)

        if len(selected_in_range) >= 10 and len(not_selected_in_range) >= 10:
            print(f"\n  重叠度区间 [{low}, {high}):")
            print(f"    选中: {len(selected_in_range)}对, 未选中: {len(not_selected_in_range)}对")

            sel_angle = sum(f['angle'] for f in selected_in_range) / len(selected_in_range)
            not_sel_angle = sum(f['angle'] for f in not_selected_in_range) / len(not_selected_in_range)
            sel_gsd = sum(f['gsd'] for f in selected_in_range) / len(selected_in_range)
            not_sel_gsd = sum(f['gsd'] for f in not_selected_in_range) / len(not_selected_in_range)

            print(f"    交会角 - 选中: {sel_angle:.2f}°, 未选中: {not_sel_angle:.2f}°")
            print(f"    GSD比 - 选中: {sel_gsd:.3f}, 未选中: {not_sel_gsd:.3f}")

    # 创新点5: 稀疏化策略分析
    print("\n" + "="*50)
    print("创新点5: 稀疏化策略分析")
    print("="*50)

    # 计算选中像对占MST的倍数
    num_nodes = len(selected_graph)
    mst_edges = num_nodes - 1
    actual_edges = len(selected_set)
    redundancy_ratio = actual_edges / mst_edges if mst_edges > 0 else 0

    print(f"  节点数: {num_nodes}")
    print(f"  MST边数: {mst_edges}")
    print(f"  实际边数: {actual_edges}")
    print(f"  冗余比: {redundancy_ratio:.2f}x MST")

    # 计算图密度
    max_edges = num_nodes * (num_nodes - 1) / 2
    density = actual_edges / max_edges if max_edges > 0 else 0
    print(f"  图密度: {density:.6f}")

    # 与随机图对比
    avg_degree = 2 * actual_edges / num_nodes if num_nodes > 0 else 0
    print(f"  平均度: {avg_degree:.2f}")

def main():
    base_dir = r"D:\codes\cpp\PairSelection\Datas"

    for i in [1, 2, 3, 4]:
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            analyze_innovation(data_dir, f"数据集{i}")

if __name__ == "__main__":
    main()
