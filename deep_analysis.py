# -*- coding: utf-8 -*-
"""
深度分析像对选择策略 - 挖掘创新点
"""

import re
import math
import os
from collections import defaultdict
import json

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

    # 也提取PhotoGroup信息
    current_photogroup = None
    photogroup_map = {}  # photo_id -> photogroup_id
    in_photogroup = False
    photogroup_id = 0

    with open(xml_path, 'r', encoding='utf-8') as f:
        for line in f:
            # PhotoGroup tracking
            if '<Photogroup>' in line or '<PhotoGroup>' in line:
                in_photogroup = True
                photogroup_id += 1
            elif '</Photogroup>' in line or '</PhotoGroup>' in line:
                in_photogroup = False

            if '<Photo>' in line:
                in_photo = True
                current_photo = {'rotation': {}, 'center': {}, 'photogroup': photogroup_id}
            elif '</Photo>' in line:
                if current_id is not None and 'x' in current_photo['center']:
                    photos[current_id] = current_photo
                    photogroup_map[current_id] = current_photo['photogroup']
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
                # 深度信息
                if '<NearDepth>' in line:
                    m = re.search(r'<NearDepth>([^<]+)</NearDepth>', line)
                    if m:
                        current_photo['near_depth'] = float(m.group(1))
                if '<MedianDepth>' in line:
                    m = re.search(r'<MedianDepth>([^<]+)</MedianDepth>', line)
                    if m:
                        current_photo['median_depth'] = float(m.group(1))
                if '<FarDepth>' in line:
                    m = re.search(r'<FarDepth>([^<]+)</FarDepth>', line)
                    if m:
                        current_photo['far_depth'] = float(m.group(1))

    return photos, photogroup_map

def get_optical_axis(photo):
    """从旋转矩阵获取光轴方向"""
    R = photo['rotation']
    if 'M_20' not in R:
        return None
    return [R['M_02'], R['M_12'], R['M_22']]

def calc_angle(v1, v2):
    """计算两个向量的夹角（度）"""
    dot = sum(a*b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a*a for a in v1))
    mag2 = math.sqrt(sum(b*b for b in v2))
    if mag1 == 0 or mag2 == 0:
        return 0
    cos_angle = max(-1, min(1, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))

def calc_distance(photo1, photo2):
    """计算空间距离"""
    c1, c2 = photo1['center'], photo2['center']
    return math.sqrt((c1['x']-c2['x'])**2 + (c1['y']-c2['y'])**2 + (c1['z']-c2['z'])**2)

def calc_convergence_angle(photo1, photo2):
    """计算交会角"""
    axis1 = get_optical_axis(photo1)
    axis2 = get_optical_axis(photo2)
    if axis1 is None or axis2 is None:
        return None
    return calc_angle(axis1, axis2)

def calc_direction_similarity(photo1, photo2):
    """计算方向相似度"""
    axis1 = get_optical_axis(photo1)
    axis2 = get_optical_axis(photo2)
    if axis1 is None or axis2 is None:
        return None
    return sum(a*b for a, b in zip(axis1, axis2))

def calc_gsd_ratio(photo1, photo2):
    """计算GSD比"""
    if 'median_depth' not in photo1 or 'median_depth' not in photo2:
        return None
    d1, d2 = photo1['median_depth'], photo2['median_depth']
    if d1 == 0 or d2 == 0:
        return None
    return max(d1, d2) / min(d1, d2)

def calc_baseline_depth_ratio(photo1, photo2):
    """计算基线/深度比 (B/D ratio) - 三角测量的重要指标"""
    distance = calc_distance(photo1, photo2)
    avg_depth = (photo1.get('median_depth', 0) + photo2.get('median_depth', 0)) / 2
    if avg_depth == 0:
        return None
    return distance / avg_depth

def parse_pairs(pairs_path):
    """解析pairs.txt"""
    pairs = []
    with open(pairs_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
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
    """从连接点提取共视关系"""
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

def percentile(values, p):
    """计算百分位数"""
    if not values:
        return 0
    values = sorted(values)
    k = (len(values) - 1) * p / 100
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return values[int(k)]
    return values[f] * (c - k) + values[c] * (k - f)

def analyze_deep(data_dir, name):
    """深度分析单个数据集"""
    print(f"\n{'='*70}")
    print(f"深度分析: {name}")
    print('='*70)

    # 找文件
    xml_file = None
    pairs_file = None
    for f in os.listdir(data_dir):
        if f.endswith('.xml'):
            xml_file = os.path.join(data_dir, f)
        if f == 'pairs.txt':
            pairs_file = os.path.join(data_dir, f)

    if not xml_file or not pairs_file:
        print("缺少必要文件")
        return None

    # 解析数据
    print("解析数据...")
    photos, photogroup_map = parse_xml_photos(xml_file)
    selected_pairs = parse_pairs(pairs_file)
    selected_set = set((min(a,b), max(a,b)) for a, b in selected_pairs)
    covisibility = parse_covisibility(xml_file)

    print(f"  照片数: {len(photos)}")
    print(f"  选中像对: {len(selected_set)}")
    print(f"  共视像对: {len(covisibility)}")

    # 构建邻接表
    neighbors = defaultdict(list)
    for a, b in selected_set:
        neighbors[a].append(b)
        neighbors[b].append(a)

    # 分析结果
    results = {
        'name': name,
        'num_photos': len(photos),
        'num_selected': len(selected_set),
        'num_covis': len(covisibility),
        'selected': [],
        'not_selected': [],
        'degree_dist': [],
        'cross_group': 0,
        'same_group': 0,
    }

    # 1. 度分布分析
    print("\n[1] 度分布分析")
    degrees = [len(neighbors[p]) for p in photos.keys() if p in neighbors]
    results['degree_dist'] = degrees
    print(f"  平均度: {sum(degrees)/len(degrees):.2f}")
    print(f"  度分布: min={min(degrees)}, max={max(degrees)}, median={sorted(degrees)[len(degrees)//2]}")

    # 度的详细分布
    degree_counts = defaultdict(int)
    for d in degrees:
        degree_counts[d] += 1
    print(f"  度频率分布(top10):")
    for d, c in sorted(degree_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    度={d}: {c}个照片 ({100*c/len(degrees):.1f}%)")

    # 2. 详细特征提取
    print("\n[2] 提取像对特征...")
    count = 0
    for pair, covis in covisibility.items():
        a, b = pair
        if a not in photos or b not in photos:
            continue

        p1, p2 = photos[a], photos[b]
        is_selected = pair in selected_set

        features = {
            'pair': pair,
            'covis': covis,
            'angle': calc_convergence_angle(p1, p2),
            'dir_sim': calc_direction_similarity(p1, p2),
            'gsd_ratio': calc_gsd_ratio(p1, p2),
            'distance': calc_distance(p1, p2),
            'bd_ratio': calc_baseline_depth_ratio(p1, p2),
            'same_group': photogroup_map.get(a, -1) == photogroup_map.get(b, -2),
        }

        if is_selected:
            results['selected'].append(features)
            if features['same_group']:
                results['same_group'] += 1
            else:
                results['cross_group'] += 1
        else:
            results['not_selected'].append(features)

        count += 1
        if count % 50000 == 0:
            print(f"  已处理 {count} 对...")

    return results

def find_decision_boundary(results):
    """分析决策边界 - 找到选中/未选中的临界条件"""
    print("\n[3] 决策边界分析")

    sel = results['selected']
    not_sel = results['not_selected']

    if not sel or not not_sel:
        return

    # 各指标的分布
    metrics = ['covis', 'angle', 'gsd_ratio', 'dir_sim', 'distance', 'bd_ratio']

    for metric in metrics:
        sel_vals = [f[metric] for f in sel if f[metric] is not None]
        not_sel_vals = [f[metric] for f in not_sel if f[metric] is not None]

        if not sel_vals or not not_sel_vals:
            continue

        # 找重叠区间
        sel_min, sel_max = min(sel_vals), max(sel_vals)
        not_sel_min, not_sel_max = min(not_sel_vals), max(not_sel_vals)

        # 计算各百分位
        sel_p25, sel_p50, sel_p75 = percentile(sel_vals, 25), percentile(sel_vals, 50), percentile(sel_vals, 75)
        not_sel_p25, not_sel_p50, not_sel_p75 = percentile(not_sel_vals, 25), percentile(not_sel_vals, 50), percentile(not_sel_vals, 75)

        print(f"\n  {metric}:")
        print(f"    选中:   P25={sel_p25:.2f}, P50={sel_p50:.2f}, P75={sel_p75:.2f}, range=[{sel_min:.2f}, {sel_max:.2f}]")
        print(f"    未选中: P25={not_sel_p25:.2f}, P50={not_sel_p50:.2f}, P75={not_sel_p75:.2f}, range=[{not_sel_min:.2f}, {not_sel_max:.2f}]")

        # 计算分离度
        if sel_p50 != not_sel_p50:
            separation = abs(sel_p50 - not_sel_p50) / max(abs(sel_p50), abs(not_sel_p50), 0.001)
            print(f"    分离度: {separation:.2f}")

def analyze_correlation(results):
    """分析因素之间的相关性"""
    print("\n[4] 因素相关性分析")

    sel = results['selected']
    if len(sel) < 100:
        print("  数据量不足")
        return

    # 提取有效数据
    data = []
    for f in sel:
        if all(f[m] is not None for m in ['covis', 'angle', 'gsd_ratio', 'dir_sim']):
            data.append(f)

    if len(data) < 100:
        return

    # 计算相关系数
    metrics = ['covis', 'angle', 'gsd_ratio', 'dir_sim']

    def correlation(x, y):
        n = len(x)
        mean_x, mean_y = sum(x)/n, sum(y)/n
        var_x = sum((xi - mean_x)**2 for xi in x)
        var_y = sum((yi - mean_y)**2 for yi in y)
        cov = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
        if var_x == 0 or var_y == 0:
            return 0
        return cov / math.sqrt(var_x * var_y)

    print("  相关系数矩阵:")
    print(f"  {'':12}", end='')
    for m in metrics:
        print(f"{m:12}", end='')
    print()

    for m1 in metrics:
        print(f"  {m1:12}", end='')
        for m2 in metrics:
            vals1 = [f[m1] for f in data]
            vals2 = [f[m2] for f in data]
            corr = correlation(vals1, vals2)
            print(f"{corr:12.3f}", end='')
        print()

def analyze_redundancy(results, photos):
    """分析冗余模式 - 每个区域有多少备用像对"""
    print("\n[5] 冗余模式分析")

    # 构建邻接关系
    neighbors = defaultdict(set)
    for f in results['selected']:
        a, b = f['pair']
        neighbors[a].add(b)
        neighbors[b].add(a)

    # 计算共同邻居数
    common_neighbor_counts = []
    for f in results['selected']:
        a, b = f['pair']
        common = neighbors[a] & neighbors[b]
        common_neighbor_counts.append(len(common))

    if common_neighbor_counts:
        avg = sum(common_neighbor_counts) / len(common_neighbor_counts)
        print(f"  像对的平均共同邻居数: {avg:.2f}")
        print(f"  分布: min={min(common_neighbor_counts)}, max={max(common_neighbor_counts)}")

        # 分布统计
        cn_counts = defaultdict(int)
        for c in common_neighbor_counts:
            cn_counts[c] += 1
        print(f"  共同邻居数分布:")
        for c in sorted(cn_counts.keys())[:10]:
            print(f"    {c}个共同邻居: {cn_counts[c]}对 ({100*cn_counts[c]/len(common_neighbor_counts):.1f}%)")

def analyze_borderline_cases(results):
    """分析边界案例 - 被选中的最差像对 vs 未被选中的最好像对"""
    print("\n[6] 边界案例分析")

    sel = results['selected']
    not_sel = results['not_selected']

    if not sel or not not_sel:
        return

    # 按共视点数排序
    sel_by_covis = sorted(sel, key=lambda x: x['covis'])
    not_sel_by_covis = sorted(not_sel, key=lambda x: -x['covis'])

    print("  被选中的最低重叠度像对 (bottom 5):")
    for f in sel_by_covis[:5]:
        print(f"    pair={f['pair']}, covis={f['covis']}, angle={f['angle']:.1f}°, gsd={f['gsd_ratio']:.2f}")

    print("\n  未被选中的最高重叠度像对 (top 5):")
    for f in not_sel_by_covis[:5]:
        print(f"    pair={f['pair']}, covis={f['covis']}, angle={f['angle']:.1f}°, gsd={f['gsd_ratio']:.2f}")

    # 分析这些边界案例的特征差异
    if len(sel_by_covis) >= 10 and len(not_sel_by_covis) >= 10:
        # 选中的最差10%
        worst_sel = sel_by_covis[:len(sel_by_covis)//10]
        # 未选中的最好10%
        best_not_sel = not_sel_by_covis[:len(not_sel_by_covis)//10]

        print("\n  边界区对比 (选中最差10% vs 未选中最好10%):")
        for metric in ['covis', 'angle', 'gsd_ratio', 'dir_sim']:
            ws_vals = [f[metric] for f in worst_sel if f[metric] is not None]
            bn_vals = [f[metric] for f in best_not_sel if f[metric] is not None]
            if ws_vals and bn_vals:
                ws_avg = sum(ws_vals) / len(ws_vals)
                bn_avg = sum(bn_vals) / len(bn_vals)
                print(f"    {metric}: 选中最差={ws_avg:.2f}, 未选中最好={bn_avg:.2f}")

def analyze_group_patterns(results):
    """分析跨组选择模式"""
    print("\n[7] 跨组选择模式")

    same = results['same_group']
    cross = results['cross_group']
    total = same + cross

    if total == 0:
        return

    print(f"  同组像对: {same} ({100*same/total:.1f}%)")
    print(f"  跨组像对: {cross} ({100*cross/total:.1f}%)")

    # 分析跨组像对的特征
    sel = results['selected']
    cross_pairs = [f for f in sel if not f['same_group']]
    same_pairs = [f for f in sel if f['same_group']]

    if cross_pairs and same_pairs:
        print("\n  跨组 vs 同组 像对特征对比:")
        for metric in ['covis', 'angle', 'gsd_ratio']:
            cross_vals = [f[metric] for f in cross_pairs if f[metric] is not None]
            same_vals = [f[metric] for f in same_pairs if f[metric] is not None]
            if cross_vals and same_vals:
                cross_avg = sum(cross_vals) / len(cross_vals)
                same_avg = sum(same_vals) / len(same_vals)
                print(f"    {metric}: 同组={same_avg:.2f}, 跨组={cross_avg:.2f}")

def analyze_selection_ratio_by_covis(results):
    """分析不同重叠度区间的选择率"""
    print("\n[8] 重叠度-选择率关系")

    sel = results['selected']
    not_sel = results['not_selected']

    # 按重叠度分桶
    buckets = [(0, 10), (10, 20), (20, 50), (50, 100), (100, 200), (200, 500), (500, 10000)]

    print(f"  {'重叠度区间':15} {'选中':8} {'未选中':8} {'选择率':10}")
    print(f"  {'-'*45}")

    for low, high in buckets:
        sel_count = sum(1 for f in sel if low <= f['covis'] < high)
        not_sel_count = sum(1 for f in not_sel if low <= f['covis'] < high)
        total = sel_count + not_sel_count
        rate = sel_count / total * 100 if total > 0 else 0
        print(f"  [{low:3d}, {high:4d})      {sel_count:8d} {not_sel_count:8d} {rate:8.1f}%")

def main():
    base_dir = r"D:\codes\cpp\PairSelection\Datas"

    all_results = []
    for i in range(1, 5):
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            result = analyze_deep(data_dir, f"数据集{i}")
            if result:
                all_results.append(result)

                # 深度分析
                find_decision_boundary(result)
                analyze_correlation(result)
                analyze_redundancy(result, None)
                analyze_borderline_cases(result)
                analyze_group_patterns(result)
                analyze_selection_ratio_by_covis(result)

    # 汇总创新点分析
    print("\n" + "="*70)
    print("创新点汇总分析")
    print("="*70)

    if all_results:
        # 合并所有数据
        all_sel = []
        all_not_sel = []
        for r in all_results:
            all_sel.extend(r['selected'])
            all_not_sel.extend(r['not_selected'])

        combined = {
            'name': '合并数据',
            'selected': all_sel,
            'not_selected': all_not_sel,
            'same_group': sum(r['same_group'] for r in all_results),
            'cross_group': sum(r['cross_group'] for r in all_results),
        }

        print("\n### 合并数据分析 ###")
        find_decision_boundary(combined)
        analyze_selection_ratio_by_covis(combined)

if __name__ == "__main__":
    main()
