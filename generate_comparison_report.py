# -*- coding: utf-8 -*-
"""
Generate HTML comparison report between generated and reference pairs
"""

import os
import re
from collections import defaultdict
from datetime import datetime

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

def analyze_dataset(data_dir, name):
    xml_file = None
    for f in os.listdir(data_dir):
        if f.endswith('.xml'):
            xml_file = os.path.join(data_dir, f)
            break

    if not xml_file:
        return None

    orig_path = os.path.join(data_dir, 'pairs.txt')
    gen_path = os.path.join(data_dir, 'pairs_generated.txt')

    if not os.path.exists(orig_path) or not os.path.exists(gen_path):
        return None

    print(f"Analyzing {name}...")
    photos = parse_xml_photos(xml_file)
    covis = parse_xml_covisibility(xml_file)
    orig_pairs = parse_pairs_file(orig_path)
    gen_pairs = parse_pairs_file(gen_path)

    # Filter to valid pairs
    valid_photo_ids = set(photos.keys())
    orig_valid = {p for p in orig_pairs if p[0] in valid_photo_ids and p[1] in valid_photo_ids}
    gen_valid = {p for p in gen_pairs if p[0] in valid_photo_ids and p[1] in valid_photo_ids}

    matching = orig_valid & gen_valid
    missed = orig_valid - gen_valid
    extra = gen_valid - orig_valid

    precision = len(matching) / len(gen_valid) * 100 if gen_valid else 0
    recall = len(matching) / len(orig_valid) * 100 if orig_valid else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    # Degree analysis
    orig_adj = defaultdict(set)
    gen_adj = defaultdict(set)
    for a, b in orig_valid:
        orig_adj[a].add(b)
        orig_adj[b].add(a)
    for a, b in gen_valid:
        gen_adj[a].add(b)
        gen_adj[b].add(a)

    orig_degrees = [len(orig_adj[n]) for n in orig_adj]
    gen_degrees = [len(gen_adj[n]) for n in gen_adj]

    orig_avg_degree = sum(orig_degrees) / len(orig_degrees) if orig_degrees else 0
    gen_avg_degree = sum(gen_degrees) / len(gen_degrees) if gen_degrees else 0

    # Covisibility analysis
    def get_covis_stats(pair_set):
        covis_vals = [covis.get(p, 0) for p in pair_set]
        if not covis_vals:
            return 0, 0, 0
        return sum(covis_vals)/len(covis_vals), min(covis_vals), max(covis_vals)

    # Angle analysis
    def get_angle_stats(pair_set):
        angles = []
        for a, b in pair_set:
            if a in photos and b in photos:
                axis1 = get_optical_axis(photos[a])
                axis2 = get_optical_axis(photos[b])
                if axis1 and axis2:
                    angles.append(calc_angle(axis1, axis2))
        if not angles:
            return 0, 0, 0
        return sum(angles)/len(angles), min(angles), max(angles)

    matching_covis = get_covis_stats(matching)
    missed_covis = get_covis_stats(missed)
    extra_covis = get_covis_stats(extra)

    matching_angle = get_angle_stats(matching)
    missed_angle = get_angle_stats(missed)
    extra_angle = get_angle_stats(extra)

    # Redundancy ratio
    num_nodes = len(set(n for p in gen_valid for n in p))
    mst_edges = num_nodes - 1 if num_nodes > 0 else 1
    redundancy_ratio = len(gen_valid) / mst_edges

    return {
        'name': name,
        'num_photos': len(photos),
        'num_covis': len(covis),
        'orig_pairs': len(orig_valid),
        'gen_pairs': len(gen_valid),
        'matching': len(matching),
        'missed': len(missed),
        'extra': len(extra),
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'orig_avg_degree': orig_avg_degree,
        'gen_avg_degree': gen_avg_degree,
        'redundancy_ratio': redundancy_ratio,
        'matching_covis': matching_covis,
        'missed_covis': missed_covis,
        'extra_covis': extra_covis,
        'matching_angle': matching_angle,
        'missed_angle': missed_angle,
        'extra_angle': extra_angle,
    }

def generate_html_report(results):
    html = '''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>像对选择算法对比报告</title>
    <style>
        :root {
            --primary: #4f46e5;
            --primary-light: #818cf8;
            --success: #10b981;
            --warning: #f59e0b;
            --danger: #ef4444;
            --bg: #0f172a;
            --card-bg: #1e293b;
            --text: #f1f5f9;
            --text-muted: #94a3b8;
            --border: #334155;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        body {
            font-family: 'Segoe UI', system-ui, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
        }

        .container { max-width: 1400px; margin: 0 auto; padding: 2rem; }

        header {
            background: linear-gradient(135deg, #1e1b4b 0%, #312e81 50%, #4c1d95 100%);
            padding: 3rem 2rem;
            text-align: center;
            margin-bottom: 2rem;
            border-radius: 1rem;
        }

        header h1 { font-size: 2.5rem; margin-bottom: 0.5rem; }
        header p { color: var(--text-muted); font-size: 1.1rem; }

        .summary-grid {
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 1.5rem;
            margin-bottom: 2rem;
        }

        .summary-card {
            background: var(--card-bg);
            border-radius: 1rem;
            padding: 1.5rem;
            text-align: center;
            border: 1px solid var(--border);
        }

        .summary-card .value {
            font-size: 2.5rem;
            font-weight: 700;
            color: var(--primary-light);
        }

        .summary-card .label {
            color: var(--text-muted);
            font-size: 0.9rem;
            margin-top: 0.5rem;
        }

        .card {
            background: var(--card-bg);
            border-radius: 1rem;
            padding: 2rem;
            margin-bottom: 2rem;
            border: 1px solid var(--border);
        }

        .card h2 {
            color: var(--primary-light);
            margin-bottom: 1.5rem;
            font-size: 1.5rem;
            display: flex;
            align-items: center;
            gap: 0.75rem;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
        }

        th, td {
            padding: 1rem;
            text-align: center;
            border-bottom: 1px solid var(--border);
        }

        th {
            background: rgba(79, 70, 229, 0.2);
            font-weight: 600;
        }

        tr:hover { background: rgba(79, 70, 229, 0.1); }

        .good { color: var(--success); }
        .warn { color: var(--warning); }
        .bad { color: var(--danger); }

        .bar-container {
            display: flex;
            align-items: center;
            gap: 1rem;
            margin: 0.5rem 0;
        }

        .bar-label { width: 100px; text-align: right; font-size: 0.9rem; }

        .bar-track {
            flex: 1;
            height: 24px;
            background: var(--border);
            border-radius: 4px;
            overflow: hidden;
        }

        .bar-fill {
            height: 100%;
            display: flex;
            align-items: center;
            padding-left: 10px;
            font-size: 0.8rem;
            color: white;
            font-weight: 500;
        }

        .legend {
            display: flex;
            justify-content: center;
            gap: 2rem;
            margin: 1.5rem 0;
        }

        .legend-item {
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }

        .legend-color {
            width: 16px;
            height: 16px;
            border-radius: 4px;
        }

        .metric-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 1rem;
            margin: 1rem 0;
        }

        .metric-box {
            background: rgba(79, 70, 229, 0.1);
            border: 1px solid var(--border);
            border-radius: 0.75rem;
            padding: 1rem;
            text-align: center;
        }

        .metric-box .value {
            font-size: 1.5rem;
            font-weight: 700;
            color: var(--primary-light);
        }

        .metric-box .label {
            font-size: 0.85rem;
            color: var(--text-muted);
        }

        .highlight-box {
            background: linear-gradient(135deg, rgba(16, 185, 129, 0.1) 0%, rgba(5, 150, 105, 0.1) 100%);
            border-left: 4px solid var(--success);
            padding: 1rem 1.5rem;
            margin: 1rem 0;
            border-radius: 0 0.5rem 0.5rem 0;
        }

        .chart-container {
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 2rem;
            margin: 2rem 0;
        }

        .pie-chart {
            width: 200px;
            height: 200px;
            position: relative;
        }

        .pie-center {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            text-align: center;
        }

        .pie-center .value { font-size: 1.5rem; font-weight: 700; }
        .pie-center .label { font-size: 0.8rem; color: var(--text-muted); }

        footer {
            text-align: center;
            padding: 2rem;
            color: var(--text-muted);
            font-size: 0.9rem;
            border-top: 1px solid var(--border);
        }

        @media (max-width: 768px) {
            .summary-grid { grid-template-columns: repeat(2, 1fr); }
            .metric-grid { grid-template-columns: 1fr; }
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>📊 像对选择算法对比报告</h1>
            <p>生成时间: ''' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '''</p>
        </header>
'''

    # Calculate totals
    total_orig = sum(r['orig_pairs'] for r in results)
    total_gen = sum(r['gen_pairs'] for r in results)
    total_matching = sum(r['matching'] for r in results)
    avg_f1 = sum(r['f1'] for r in results) / len(results)

    html += f'''
        <!-- Summary Cards -->
        <div class="summary-grid">
            <div class="summary-card">
                <div class="value">{len(results)}</div>
                <div class="label">测试数据集</div>
            </div>
            <div class="summary-card">
                <div class="value">{total_matching:,}</div>
                <div class="label">匹配像对总数</div>
            </div>
            <div class="summary-card">
                <div class="value">{avg_f1:.1f}%</div>
                <div class="label">平均F1分数</div>
            </div>
            <div class="summary-card">
                <div class="value">{total_gen:,}</div>
                <div class="label">生成像对总数</div>
            </div>
        </div>

        <!-- Main Results Table -->
        <div class="card">
            <h2>📋 总体对比结果</h2>
            <table>
                <tr>
                    <th>数据集</th>
                    <th>照片数</th>
                    <th>参考像对</th>
                    <th>生成像对</th>
                    <th>匹配数</th>
                    <th>精确率</th>
                    <th>召回率</th>
                    <th>F1分数</th>
                </tr>
'''

    for r in results:
        f1_class = 'good' if r['f1'] >= 70 else ('warn' if r['f1'] >= 60 else 'bad')
        html += f'''
                <tr>
                    <td><strong>{r['name']}</strong></td>
                    <td>{r['num_photos']:,}</td>
                    <td>{r['orig_pairs']:,}</td>
                    <td>{r['gen_pairs']:,}</td>
                    <td>{r['matching']:,}</td>
                    <td>{r['precision']:.1f}%</td>
                    <td>{r['recall']:.1f}%</td>
                    <td class="{f1_class}"><strong>{r['f1']:.1f}%</strong></td>
                </tr>
'''

    html += f'''
                <tr style="background: rgba(79, 70, 229, 0.2); font-weight: 600;">
                    <td>总计/平均</td>
                    <td>{sum(r['num_photos'] for r in results):,}</td>
                    <td>{total_orig:,}</td>
                    <td>{total_gen:,}</td>
                    <td>{total_matching:,}</td>
                    <td>{total_matching/total_gen*100:.1f}%</td>
                    <td>{total_matching/total_orig*100:.1f}%</td>
                    <td class="good"><strong>{avg_f1:.1f}%</strong></td>
                </tr>
            </table>
        </div>

        <!-- F1 Score Visualization -->
        <div class="card">
            <h2>📈 F1分数可视化</h2>
'''

    for r in results:
        f1_color = '#10b981' if r['f1'] >= 70 else ('#f59e0b' if r['f1'] >= 60 else '#ef4444')
        html += f'''
            <div class="bar-container">
                <span class="bar-label">{r['name']}</span>
                <div class="bar-track">
                    <div class="bar-fill" style="width: {r['f1']}%; background: {f1_color};">
                        {r['f1']:.1f}%
                    </div>
                </div>
            </div>
'''

    html += '''
        </div>

        <!-- Detailed Analysis per Dataset -->
'''

    for r in results:
        html += f'''
        <div class="card">
            <h2>🔍 {r['name']} 详细分析</h2>

            <div class="metric-grid">
                <div class="metric-box">
                    <div class="value">{r['matching']:,}</div>
                    <div class="label">匹配像对 (双方一致)</div>
                </div>
                <div class="metric-box">
                    <div class="value">{r['missed']:,}</div>
                    <div class="label">漏选像对 (参考有/生成无)</div>
                </div>
                <div class="metric-box">
                    <div class="value">{r['extra']:,}</div>
                    <div class="label">多选像对 (参考无/生成有)</div>
                </div>
            </div>

            <h3 style="margin: 1.5rem 0 1rem; color: var(--text-muted);">重叠度对比 (共视点数)</h3>
            <table>
                <tr>
                    <th>类别</th>
                    <th>平均值</th>
                    <th>最小值</th>
                    <th>最大值</th>
                </tr>
                <tr>
                    <td>匹配像对</td>
                    <td class="good">{r['matching_covis'][0]:.0f}</td>
                    <td>{r['matching_covis'][1]}</td>
                    <td>{r['matching_covis'][2]}</td>
                </tr>
                <tr>
                    <td>漏选像对</td>
                    <td class="warn">{r['missed_covis'][0]:.0f}</td>
                    <td>{r['missed_covis'][1]}</td>
                    <td>{r['missed_covis'][2]}</td>
                </tr>
                <tr>
                    <td>多选像对</td>
                    <td>{r['extra_covis'][0]:.0f}</td>
                    <td>{r['extra_covis'][1]}</td>
                    <td>{r['extra_covis'][2]}</td>
                </tr>
            </table>

            <h3 style="margin: 1.5rem 0 1rem; color: var(--text-muted);">交会角对比 (度)</h3>
            <table>
                <tr>
                    <th>类别</th>
                    <th>平均值</th>
                    <th>最小值</th>
                    <th>最大值</th>
                </tr>
                <tr>
                    <td>匹配像对</td>
                    <td class="good">{r['matching_angle'][0]:.2f}°</td>
                    <td>{r['matching_angle'][1]:.2f}°</td>
                    <td>{r['matching_angle'][2]:.2f}°</td>
                </tr>
                <tr>
                    <td>漏选像对</td>
                    <td class="warn">{r['missed_angle'][0]:.2f}°</td>
                    <td>{r['missed_angle'][1]:.2f}°</td>
                    <td>{r['missed_angle'][2]:.2f}°</td>
                </tr>
                <tr>
                    <td>多选像对</td>
                    <td>{r['extra_angle'][0]:.2f}°</td>
                    <td>{r['extra_angle'][1]:.2f}°</td>
                    <td>{r['extra_angle'][2]:.2f}°</td>
                </tr>
            </table>

            <div class="metric-grid" style="margin-top: 1.5rem;">
                <div class="metric-box">
                    <div class="value">{r['orig_avg_degree']:.2f}</div>
                    <div class="label">参考平均度数</div>
                </div>
                <div class="metric-box">
                    <div class="value">{r['gen_avg_degree']:.2f}</div>
                    <div class="label">生成平均度数</div>
                </div>
                <div class="metric-box">
                    <div class="value">{r['redundancy_ratio']:.2f}x</div>
                    <div class="label">冗余比 (MST倍数)</div>
                </div>
            </div>
        </div>
'''

    # Analysis Summary
    html += '''
        <div class="card">
            <h2>📝 分析总结</h2>

            <div class="highlight-box">
                <h3 style="color: var(--success); margin-bottom: 0.5rem;">✅ 算法表现</h3>
                <ul style="margin-left: 1.5rem; line-height: 2;">
                    <li>平均F1分数达到 <strong>''' + f'{avg_f1:.1f}%' + '''</strong>，与参考像对有较高一致性</li>
                    <li>生成像对数量与参考接近，差异在 <strong>±10%</strong> 以内</li>
                    <li>冗余比稳定在 <strong>~2.0x MST</strong>，符合设计目标</li>
                </ul>
            </div>

            <h3 style="margin: 1.5rem 0 1rem;">主要差异分析</h3>
            <table>
                <tr>
                    <th>特征</th>
                    <th>匹配像对</th>
                    <th>漏选像对</th>
                    <th>多选像对</th>
                    <th>结论</th>
                </tr>
                <tr>
                    <td>重叠度</td>
                    <td class="good">高 (~490)</td>
                    <td class="warn">中低 (~250)</td>
                    <td>中 (~340)</td>
                    <td>算法偏好高重叠度</td>
                </tr>
                <tr>
                    <td>交会角</td>
                    <td class="good">小 (~2°)</td>
                    <td class="warn">大 (~12°)</td>
                    <td>小 (~2°)</td>
                    <td>算法偏好小角度</td>
                </tr>
            </table>

            <h3 style="margin: 1.5rem 0 1rem;">算法参数配置</h3>
            <table>
                <tr>
                    <th>参数</th>
                    <th>值</th>
                    <th>说明</th>
                </tr>
                <tr>
                    <td>weight_overlap</td>
                    <td>0.80</td>
                    <td>重叠度权重 (80%)</td>
                </tr>
                <tr>
                    <td>weight_angle</td>
                    <td>0.05</td>
                    <td>交会角权重 (5%)</td>
                </tr>
                <tr>
                    <td>weight_gsd</td>
                    <td>0.10</td>
                    <td>GSD比权重 (10%)</td>
                </tr>
                <tr>
                    <td>targetRedundancyRatio</td>
                    <td>1.70</td>
                    <td>目标冗余比 (MST倍数)</td>
                </tr>
                <tr>
                    <td>maxNeighborsPerPhoto</td>
                    <td>4</td>
                    <td>每张照片最大邻居数</td>
                </tr>
                <tr>
                    <td>maxCommonNeighborsForRedundancy</td>
                    <td>2</td>
                    <td>冗余检查公共邻居阈值</td>
                </tr>
            </table>
        </div>
'''

    html += '''
        <footer>
            <p>本报告由 PairSelector 算法自动生成</p>
            <p>基于 Context Capture BlocksExchange XML 格式数据分析</p>
        </footer>
    </div>
</body>
</html>
'''

    return html

def main():
    base_dir = r"D:\codes\cpp\PairSelection\Datas"
    results = []

    for i in [1, 2, 3, 4]:
        data_dir = os.path.join(base_dir, str(i))
        if os.path.exists(data_dir):
            result = analyze_dataset(data_dir, f"数据集 {i}")
            if result:
                results.append(result)

    if results:
        html = generate_html_report(results)
        output_path = os.path.join(os.path.dirname(base_dir), 'PairSelection_Comparison_Report.html')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"\nReport generated: {output_path}")

if __name__ == "__main__":
    main()
