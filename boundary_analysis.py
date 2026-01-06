# -*- coding: utf-8 -*-
"""
深入分析边界案例和决策逻辑
探索为什么高重叠度像对没被选中
"""

import re
import math
import os
from collections import defaultdict

def parse_xml_photos(xml_path, max_photos=5000):
    """从XML提取照片信息"""
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
                current_photo = {}
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

def parse_covisibility(xml_path, max_ties=50000):
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

def analyze_why_rejected(data_dir, name):
    """分析为什么高重叠度像对被拒绝"""
    print(f"\n{'='*70}")
    print(f"边界案例深度分析: {name}")
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

    # 构建每张照片的邻居列表和选中邻居列表
    all_neighbors = defaultdict(list)  # photo -> [(neighbor, covis, angle, gsd)]
    selected_neighbors = defaultdict(list)

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

    # 分析每张照片的选择模式
    print("\n[1] 每张照片的选择模式分析")

    # 统计: 每张照片选了多少邻居，这些邻居在所有候选中排名如何
    rank_stats = []
    top_k_selection = []

    for photo_id, neighbors in all_neighbors.items():
        if len(neighbors) < 5:
            continue

        # 按重叠度排序
        sorted_by_covis = sorted(neighbors, key=lambda x: -x['covis'])
        selected_count = sum(1 for n in neighbors if n['selected'])

        # 计算选中邻居的平均排名
        selected_ranks = []
        for rank, n in enumerate(sorted_by_covis):
            if n['selected']:
                selected_ranks.append(rank + 1)

        if selected_ranks:
            avg_rank = sum(selected_ranks) / len(selected_ranks)
            rank_stats.append({
                'photo': photo_id,
                'total_neighbors': len(neighbors),
                'selected_count': selected_count,
                'avg_rank': avg_rank,
                'max_rank': max(selected_ranks),
                'ranks': selected_ranks
            })

            # 检查是否严格选择Top-K
            is_strict_topk = all(r <= selected_count + 1 for r in selected_ranks)
            top_k_selection.append(is_strict_topk)

    if rank_stats:
        avg_avg_rank = sum(r['avg_rank'] for r in rank_stats) / len(rank_stats)
        avg_max_rank = sum(r['max_rank'] for r in rank_stats) / len(rank_stats)
        strict_topk_rate = sum(top_k_selection) / len(top_k_selection) * 100

        print(f"  分析了 {len(rank_stats)} 张照片")
        print(f"  选中邻居的平均排名: {avg_avg_rank:.2f}")
        print(f"  选中邻居的最大排名平均值: {avg_max_rank:.2f}")
        print(f"  严格Top-K选择率: {strict_topk_rate:.1f}%")

    # 分析"反常"案例：高重叠但未选中
    print("\n[2] 高重叠度但未被选中的案例分析")

    rejected_high_covis = []
    for pair, covis in covisibility.items():
        if pair in selected_set:
            continue
        a, b = pair
        if a not in photos or b not in photos:
            continue
        if covis < 200:  # 只看高重叠度的
            continue

        p1, p2 = photos[a], photos[b]
        angle = calc_convergence_angle(p1, p2)
        gsd = calc_gsd_ratio(p1, p2)

        # 检查a和b各自选了谁
        a_selected = [n for n in all_neighbors[a] if n['selected']]
        b_selected = [n for n in all_neighbors[b] if n['selected']]

        # 检查这对在各自的排名
        a_neighbors_sorted = sorted(all_neighbors[a], key=lambda x: -x['covis'])
        b_neighbors_sorted = sorted(all_neighbors[b], key=lambda x: -x['covis'])

        a_rank = next((i+1 for i, n in enumerate(a_neighbors_sorted) if n['id'] == b), -1)
        b_rank = next((i+1 for i, n in enumerate(b_neighbors_sorted) if n['id'] == a), -1)

        rejected_high_covis.append({
            'pair': pair,
            'covis': covis,
            'angle': angle,
            'gsd': gsd,
            'a_rank': a_rank,
            'b_rank': b_rank,
            'a_selected_count': len(a_selected),
            'b_selected_count': len(b_selected),
            'a_total': len(all_neighbors[a]),
            'b_total': len(all_neighbors[b]),
        })

    rejected_high_covis.sort(key=lambda x: -x['covis'])

    print(f"  高重叠度(>=200)但未选中的像对数: {len(rejected_high_covis)}")
    print("\n  典型案例 (按重叠度排序前10):")
    for case in rejected_high_covis[:10]:
        print(f"    pair={case['pair']}, covis={case['covis']}")
        print(f"      angle={case['angle']:.1f}°, gsd={case['gsd']:.2f}")
        print(f"      在照片{case['pair'][0]}中排名: {case['a_rank']}/{case['a_total']}, 该照片已选{case['a_selected_count']}个邻居")
        print(f"      在照片{case['pair'][1]}中排名: {case['b_rank']}/{case['b_total']}, 该照片已选{case['b_selected_count']}个邻居")

    # 分析这些被拒绝案例的共同点
    if rejected_high_covis:
        avg_a_rank = sum(c['a_rank'] for c in rejected_high_covis) / len(rejected_high_covis)
        avg_b_rank = sum(c['b_rank'] for c in rejected_high_covis) / len(rejected_high_covis)
        print(f"\n  被拒绝高重叠像对的平均排名: A方={avg_a_rank:.1f}, B方={avg_b_rank:.1f}")

        # 排名分布
        rank_dist = defaultdict(int)
        for c in rejected_high_covis:
            max_rank = max(c['a_rank'], c['b_rank'])
            if max_rank <= 5:
                rank_dist['<=5'] += 1
            elif max_rank <= 10:
                rank_dist['6-10'] += 1
            elif max_rank <= 20:
                rank_dist['11-20'] += 1
            else:
                rank_dist['>20'] += 1

        print(f"  最大排名分布:")
        for k, v in sorted(rank_dist.items()):
            print(f"    {k}: {v}对 ({100*v/len(rejected_high_covis):.1f}%)")

    # 分析选择是否考虑"覆盖均匀性"
    print("\n[3] 覆盖均匀性分析")

    # 检查是否存在"区域配额"模式
    # 如果两张照片已经通过其他路径连接，可能不需要直接连接

    # 计算图的路径
    def bfs_distance(graph, start, max_dist=3):
        """BFS计算从start到其他节点的距离"""
        dist = {start: 0}
        queue = [start]
        while queue:
            node = queue.pop(0)
            if dist[node] >= max_dist:
                continue
            for neighbor in graph.get(node, []):
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    queue.append(neighbor)
        return dist

    # 构建选中像对的图
    selected_graph = defaultdict(list)
    for pair in selected_set:
        a, b = pair
        selected_graph[a].append(b)
        selected_graph[b].append(a)

    # 检查被拒绝的高重叠像对，是否在图中已经有短路径
    if rejected_high_covis:
        path_lengths = []
        for case in rejected_high_covis[:50]:  # 只检查前50个
            a, b = case['pair']
            distances = bfs_distance(selected_graph, a, max_dist=4)
            if b in distances:
                path_lengths.append(distances[b])
            else:
                path_lengths.append(99)  # 不可达

        connected_2hop = sum(1 for p in path_lengths if p == 2)
        connected_3hop = sum(1 for p in path_lengths if p == 3)
        disconnected = sum(1 for p in path_lengths if p > 3)

        print(f"  被拒绝的高重叠像对在选中图中的连通性:")
        print(f"    2-hop连通: {connected_2hop} ({100*connected_2hop/len(path_lengths):.1f}%)")
        print(f"    3-hop连通: {connected_3hop} ({100*connected_3hop/len(path_lengths):.1f}%)")
        print(f"    >3-hop或不连通: {disconnected} ({100*disconnected/len(path_lengths):.1f}%)")

    # 分析选择是否有"多样性"考虑
    print("\n[4] 选择多样性分析")

    # 检查每张照片选中的邻居是否来自不同方向
    diversity_scores = []
    for photo_id, neighbors in all_neighbors.items():
        selected = [n for n in neighbors if n['selected']]
        if len(selected) < 2:
            continue

        # 计算选中邻居之间的交会角多样性
        angles_between_selected = []
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                n1_id, n2_id = selected[i]['id'], selected[j]['id']
                if n1_id in photos and n2_id in photos:
                    angle = calc_convergence_angle(photos[n1_id], photos[n2_id])
                    if angle is not None:
                        angles_between_selected.append(angle)

        if angles_between_selected:
            avg_diversity = sum(angles_between_selected) / len(angles_between_selected)
            diversity_scores.append(avg_diversity)

    if diversity_scores:
        print(f"  选中邻居之间的平均交会角: {sum(diversity_scores)/len(diversity_scores):.2f}°")
        print(f"  范围: [{min(diversity_scores):.1f}°, {max(diversity_scores):.1f}°]")

def main():
    base_dir = r"D:\codes\cpp\PairSelection\Datas"

    for i in [1, 4]:  # 分析小数据集和大数据集
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            analyze_why_rejected(data_dir, f"数据集{i}")

if __name__ == "__main__":
    main()
