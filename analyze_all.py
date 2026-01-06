# -*- coding: utf-8 -*-
"""
分析4个数据集的像对选择策略共同点
考虑因素：重叠度、交会角、分辨率、方向向量
"""

import re
import math
import os
from collections import defaultdict

def parse_xml_photos(xml_path, max_photos=5000):
    """从XML提取照片信息：位姿、旋转矩阵、深度"""
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

    return photos

def get_optical_axis(photo):
    """从旋转矩阵获取光轴方向（相机Z轴在世界坐标系中的方向）"""
    R = photo['rotation']
    if 'M_20' not in R:
        return None
    # 光轴是旋转矩阵的第三列（或第三行，取决于约定）
    # Context Capture使用的是R*X_cam = X_world，所以光轴是第三列
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

def calc_convergence_angle(photo1, photo2):
    """计算交会角：两个相机光轴指向对方中心形成的角度"""
    c1 = photo1['center']
    c2 = photo2['center']

    # 基线向量（从相机1指向相机2）
    baseline = [c2['x'] - c1['x'], c2['y'] - c1['y'], c2['z'] - c1['z']]
    baseline_len = math.sqrt(sum(b*b for b in baseline))
    if baseline_len == 0:
        return 0
    baseline = [b/baseline_len for b in baseline]

    # 获取两个相机的光轴
    axis1 = get_optical_axis(photo1)
    axis2 = get_optical_axis(photo2)

    if axis1 is None or axis2 is None:
        return None

    # 交会角 = 两个光轴方向的夹角
    angle = calc_angle(axis1, axis2)
    return angle

def calc_direction_similarity(photo1, photo2):
    """计算方向向量相似度（光轴方向的余弦相似度）"""
    axis1 = get_optical_axis(photo1)
    axis2 = get_optical_axis(photo2)

    if axis1 is None or axis2 is None:
        return None

    dot = sum(a*b for a, b in zip(axis1, axis2))
    return dot  # 余弦值，1表示同向，-1表示反向，0表示正交

def calc_gsd_ratio(photo1, photo2):
    """计算分辨率比（通过中值深度比来近似GSD比）"""
    if 'median_depth' not in photo1 or 'median_depth' not in photo2:
        return None
    d1 = photo1['median_depth']
    d2 = photo2['median_depth']
    if d1 == 0 or d2 == 0:
        return None
    ratio = max(d1, d2) / min(d1, d2)
    return ratio

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

def parse_covisibility(xml_path, max_ties=20000):
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

def analyze_dataset(data_dir, name):
    """分析单个数据集"""
    print(f"\n{'='*60}")
    print(f"分析数据集: {name}")
    print('='*60)

    # 找到XML和pairs文件
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
    print("解析照片信息...")
    photos = parse_xml_photos(xml_file)
    print(f"  照片数: {len(photos)}")

    print("解析选中的像对...")
    selected_pairs = parse_pairs(pairs_file)
    selected_set = set((min(a,b), max(a,b)) for a, b in selected_pairs)
    print(f"  选中像对数: {len(selected_set)}")

    print("解析共视关系...")
    covisibility = parse_covisibility(xml_file)
    print(f"  共视像对数: {len(covisibility)}")

    # 分析各因素
    results = {
        'name': name,
        'num_photos': len(photos),
        'num_selected': len(selected_set),
        'num_covis': len(covisibility),
        'convergence_angles': {'selected': [], 'not_selected': []},
        'direction_similarity': {'selected': [], 'not_selected': []},
        'gsd_ratio': {'selected': [], 'not_selected': []},
        'covis_count': {'selected': [], 'not_selected': []},
    }

    print("计算各因素...")
    count = 0
    for pair, covis in covisibility.items():
        a, b = pair
        if a not in photos or b not in photos:
            continue

        p1, p2 = photos[a], photos[b]
        is_selected = pair in selected_set

        # 交会角
        angle = calc_convergence_angle(p1, p2)
        if angle is not None:
            if is_selected:
                results['convergence_angles']['selected'].append(angle)
            else:
                results['convergence_angles']['not_selected'].append(angle)

        # 方向相似度
        sim = calc_direction_similarity(p1, p2)
        if sim is not None:
            if is_selected:
                results['direction_similarity']['selected'].append(sim)
            else:
                results['direction_similarity']['not_selected'].append(sim)

        # GSD比
        gsd = calc_gsd_ratio(p1, p2)
        if gsd is not None:
            if is_selected:
                results['gsd_ratio']['selected'].append(gsd)
            else:
                results['gsd_ratio']['not_selected'].append(gsd)

        # 共视点数
        if is_selected:
            results['covis_count']['selected'].append(covis)
        else:
            results['covis_count']['not_selected'].append(covis)

        count += 1
        if count % 50000 == 0:
            print(f"  已处理 {count} 对...")

    return results

def print_stats(values, name):
    """打印统计信息"""
    if not values:
        print(f"  {name}: 无数据")
        return
    values = sorted(values)
    n = len(values)
    print(f"  {name}:")
    print(f"    数量: {n}")
    print(f"    平均: {sum(values)/n:.2f}")
    print(f"    中位数: {values[n//2]:.2f}")
    print(f"    最小: {min(values):.2f}")
    print(f"    最大: {max(values):.2f}")
    print(f"    25%: {values[n//4]:.2f}")
    print(f"    75%: {values[3*n//4]:.2f}")

def main():
    base_dir = r"D:\codes\cpp\PairSelection\Datas"

    all_results = []
    for i in range(1, 5):
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            result = analyze_dataset(data_dir, f"数据集{i}")
            if result:
                all_results.append(result)

    # 打印各数据集的详细统计
    print("\n" + "="*80)
    print("各数据集统计对比")
    print("="*80)

    for r in all_results:
        print(f"\n### {r['name']} ###")
        print(f"照片数: {r['num_photos']}, 选中像对: {r['num_selected']}, 共视像对: {r['num_covis']}")

        print("\n[交会角(度)]")
        print_stats(r['convergence_angles']['selected'], "选中")
        print_stats(r['convergence_angles']['not_selected'], "未选中")

        print("\n[方向相似度(余弦值)]")
        print_stats(r['direction_similarity']['selected'], "选中")
        print_stats(r['direction_similarity']['not_selected'], "未选中")

        print("\n[GSD比(分辨率比)]")
        print_stats(r['gsd_ratio']['selected'], "选中")
        print_stats(r['gsd_ratio']['not_selected'], "未选中")

        print("\n[共视点数(重叠度)]")
        print_stats(r['covis_count']['selected'], "选中")
        print_stats(r['covis_count']['not_selected'], "未选中")

    # 汇总共同点
    print("\n" + "="*80)
    print("共同点分析汇总")
    print("="*80)

    # 合并所有数据集的统计
    all_selected = {'angle': [], 'dir_sim': [], 'gsd': [], 'covis': []}
    all_not_selected = {'angle': [], 'dir_sim': [], 'gsd': [], 'covis': []}

    for r in all_results:
        all_selected['angle'].extend(r['convergence_angles']['selected'])
        all_selected['dir_sim'].extend(r['direction_similarity']['selected'])
        all_selected['gsd'].extend(r['gsd_ratio']['selected'])
        all_selected['covis'].extend(r['covis_count']['selected'])

        all_not_selected['angle'].extend(r['convergence_angles']['not_selected'])
        all_not_selected['dir_sim'].extend(r['direction_similarity']['not_selected'])
        all_not_selected['gsd'].extend(r['gsd_ratio']['not_selected'])
        all_not_selected['covis'].extend(r['covis_count']['not_selected'])

    print("\n汇总统计(所有数据集合并):")

    for metric, name in [('angle', '交会角'), ('dir_sim', '方向相似度'), ('gsd', 'GSD比'), ('covis', '共视点数')]:
        sel = all_selected[metric]
        not_sel = all_not_selected[metric]
        if sel and not_sel:
            sel_avg = sum(sel) / len(sel)
            not_sel_avg = sum(not_sel) / len(not_sel)
            print(f"\n{name}:")
            print(f"  选中平均: {sel_avg:.2f}")
            print(f"  未选中平均: {not_sel_avg:.2f}")
            print(f"  差异: {sel_avg - not_sel_avg:.2f} ({(sel_avg/not_sel_avg - 1)*100:.1f}%)")

if __name__ == "__main__":
    main()
