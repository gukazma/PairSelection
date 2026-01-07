# -*- coding: utf-8 -*-
"""
深度分析像对选择策略
根据pairs.txt三段式结构：
1. 密集匹配对 (Dense Matching Pairs)
2. Refine对 (Refinement Pairs)
3. 密集匹配三角形 (Dense Matching Triplets)
"""

import os
import re
from collections import defaultdict

def parse_pairs_file_sections(path):
    """解析pairs.txt，分离三个区段"""
    sections = []
    current_section = []

    with open(path, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1

        # 跳过注释
        if not line or line.startswith('#'):
            continue

        parts = line.split()

        # 单个数字 = 区段计数器
        if len(parts) == 1:
            if current_section:
                sections.append(current_section)
            current_section = []
            continue

        # 2个或3个数字 = 数据行
        if len(parts) >= 2:
            try:
                ids = [int(p) for p in parts]
                current_section.append(ids)
            except:
                pass

    if current_section:
        sections.append(current_section)

    return sections

def extract_pairs_from_section(section):
    """从区段提取所有像对"""
    pairs = set()
    for ids in section:
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = min(ids[i], ids[j]), max(ids[i], ids[j])
                pairs.add((a, b))
    return pairs

def parse_generated_pairs(path):
    """解析生成的pairs文件"""
    pairs = set()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 2:
                try:
                    a, b = int(parts[0]), int(parts[1])
                    pairs.add((min(a,b), max(a,b)))
                except:
                    pass
    return pairs

def analyze_dataset(data_dir, name):
    """分析单个数据集"""
    pairs_path = os.path.join(data_dir, 'pairs.txt')
    gen_path = os.path.join(data_dir, 'pairs_generated.txt')

    if not os.path.exists(pairs_path) or not os.path.exists(gen_path):
        return None

    # 解析原始文件的三个区段
    sections = parse_pairs_file_sections(pairs_path)

    print(f"\n{'='*60}")
    print(f"数据集: {name}")
    print(f"{'='*60}")

    if len(sections) >= 1:
        dense_pairs = extract_pairs_from_section(sections[0])
        print(f"\n第一区段 (密集匹配对):")
        print(f"  - 原始行数: {len(sections[0])}")
        print(f"  - 唯一像对: {len(dense_pairs)}")
    else:
        dense_pairs = set()

    if len(sections) >= 2:
        refine_pairs = extract_pairs_from_section(sections[1])
        print(f"\n第二区段 (Refine对):")
        print(f"  - 原始行数: {len(sections[1])}")
        print(f"  - 唯一像对: {len(refine_pairs)}")
    else:
        refine_pairs = set()

    if len(sections) >= 3:
        triplet_pairs = extract_pairs_from_section(sections[2])
        triplet_count = len(sections[2])
        print(f"\n第三区段 (密集匹配三角形):")
        print(f"  - 三角形数: {triplet_count}")
        print(f"  - 展开像对: {len(triplet_pairs)}")
    else:
        triplet_pairs = set()
        triplet_count = 0

    # 计算各种组合
    dense_all = dense_pairs | triplet_pairs  # 所有密集匹配像对
    all_pairs = dense_pairs | refine_pairs | triplet_pairs  # 所有像对

    print(f"\n组合分析:")
    print(f"  - 密集匹配总像对 (区段1+区段3): {len(dense_all)}")
    print(f"  - 全部像对 (区段1+2+3): {len(all_pairs)}")

    # 区段重叠分析
    overlap_1_2 = dense_pairs & refine_pairs
    overlap_1_3 = dense_pairs & triplet_pairs
    overlap_2_3 = refine_pairs & triplet_pairs

    print(f"\n区段重叠:")
    print(f"  - 区段1 ∩ 区段2: {len(overlap_1_2)} 对")
    print(f"  - 区段1 ∩ 区段3: {len(overlap_1_3)} 对")
    print(f"  - 区段2 ∩ 区段3: {len(overlap_2_3)} 对")

    # 加载生成的像对
    gen_pairs = parse_generated_pairs(gen_path)
    print(f"\n生成像对: {len(gen_pairs)}")

    # 与各区段比较
    print(f"\n与生成像对的匹配情况:")

    def calc_metrics(ref, gen):
        match = ref & gen
        prec = len(match)/len(gen)*100 if gen else 0
        rec = len(match)/len(ref)*100 if ref else 0
        f1 = 2*prec*rec/(prec+rec) if (prec+rec) > 0 else 0
        return len(match), prec, rec, f1

    # 与第一区段比较
    m, p, r, f = calc_metrics(dense_pairs, gen_pairs)
    print(f"  vs 区段1 (密集匹配对):     匹配={m:4d}, P={p:5.1f}%, R={r:5.1f}%, F1={f:5.1f}%")

    # 与第二区段比较
    m, p, r, f = calc_metrics(refine_pairs, gen_pairs)
    print(f"  vs 区段2 (Refine对):       匹配={m:4d}, P={p:5.1f}%, R={r:5.1f}%, F1={f:5.1f}%")

    # 与第三区段展开比较
    m, p, r, f = calc_metrics(triplet_pairs, gen_pairs)
    print(f"  vs 区段3 (三角形展开):     匹配={m:4d}, P={p:5.1f}%, R={r:5.1f}%, F1={f:5.1f}%")

    # 与密集匹配总体比较
    m, p, r, f = calc_metrics(dense_all, gen_pairs)
    print(f"  vs 密集匹配总体 (1+3):     匹配={m:4d}, P={p:5.1f}%, R={r:5.1f}%, F1={f:5.1f}%")

    # 与全部比较
    m, p, r, f = calc_metrics(all_pairs, gen_pairs)
    print(f"  vs 全部像对 (1+2+3):       匹配={m:4d}, P={p:5.1f}%, R={r:5.1f}%, F1={f:5.1f}%")

    # 分析生成像对的覆盖情况
    only_in_1 = gen_pairs & dense_pairs - refine_pairs - triplet_pairs
    only_in_2 = gen_pairs & refine_pairs - dense_pairs - triplet_pairs
    only_in_3 = gen_pairs & triplet_pairs - dense_pairs - refine_pairs
    in_multiple = gen_pairs & (dense_pairs & refine_pairs | dense_pairs & triplet_pairs | refine_pairs & triplet_pairs)
    not_in_ref = gen_pairs - all_pairs

    print(f"\n生成像对来源分布:")
    print(f"  - 仅在区段1: {len(only_in_1)}")
    print(f"  - 仅在区段2: {len(only_in_2)}")
    print(f"  - 仅在区段3: {len(only_in_3)}")
    print(f"  - 多区段重叠: {len(in_multiple)}")
    print(f"  - 不在参考中: {len(not_in_ref)}")

    return {
        'name': name,
        'dense_pairs': len(dense_pairs),
        'refine_pairs': len(refine_pairs),
        'triplet_pairs': len(triplet_pairs),
        'dense_all': len(dense_all),
        'all_pairs': len(all_pairs),
        'gen_pairs': len(gen_pairs),
        'sections': sections
    }

def main():
    base_dir = "D:/codes/cpp/PairSelection/Datas"

    print("="*60)
    print("像对选择策略深度分析")
    print("="*60)

    results = []
    for i in range(1, 5):
        data_dir = os.path.join(base_dir, str(i))
        result = analyze_dataset(data_dir, f"数据集 {i}")
        if result:
            results.append(result)

    # 总结
    print(f"\n{'='*60}")
    print("总结")
    print(f"{'='*60}")

    print(f"\n各数据集像对数量:")
    print(f"{'数据集':<10} {'密集对':<8} {'Refine对':<10} {'三角形对':<10} {'密集总计':<10} {'全部':<8} {'生成':<8}")
    print("-"*70)
    for r in results:
        print(f"{r['name']:<10} {r['dense_pairs']:<8} {r['refine_pairs']:<10} {r['triplet_pairs']:<10} {r['dense_all']:<10} {r['all_pairs']:<8} {r['gen_pairs']:<8}")

if __name__ == "__main__":
    main()
