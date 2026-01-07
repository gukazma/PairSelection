# 像对选择算法研发方案 v2.0

## 1. 问题分析

### 1.1 参考输出结构

参考 `pairs.txt` 包含三个独立区段：

```
#Graph of Views
#Version 2
[count1]          ← 第一区段: 密集匹配对数量
id1 id2           ← 密集匹配像对 (Dense Matching Pairs)
...
[count2]          ← 第二区段: Refine对数量
id1 id2           ← Refine像对 (Refinement Pairs)
...
[count3]          ← 第三区段: 三角形数量
id1 id2 id3       ← 三角形约束 (Triplet Constraints)
...
```

### 1.2 关键发现

| 特性 | 观察结果 |
|------|----------|
| 区段1 ⊂ 区段3 | 密集匹配对100%包含在三角形展开边中 |
| 区段1 ∩ 区段2 | ~85-92%重叠 |
| 三角形核心 | 三角形是整个策略的核心，边从三角形导出 |
| 冗余比 | 三角形展开后约1.8-2.0x MST |

### 1.3 各区段统计 (参考数据)

| 数据集 | 照片数 | 区段1 (密集) | 区段2 (Refine) | 区段3 (三角形) | 三角形展开 |
|--------|--------|-------------|---------------|---------------|-----------|
| 1 | 175 | 198 | 198 | 111 | 309 |
| 2 | 340 | 393 | 400 | 223 | 616 |
| 3 | 300 | 359 | 361 | 205 | 561 |
| 4 | 3425 | 4226 | 4341 | 2353 | 6479 |

### 1.4 当前算法问题

- 只输出单一像对列表，没有三段结构
- 没有显式三角形约束输出
- 无法区分密集匹配和Refine用途
- F1约73%，但召回率不高

## 2. 新算法设计

### 2.1 核心思路

**三角形优先策略**: 先选择高质量三角形，再从三角形边中导出密集匹配对。

```
候选像对 → 构建候选三角形 → 选择优质三角形 → 导出密集对 → 生成Refine对
```

### 2.2 算法流程

```
Phase 1: 数据准备
├── 解析XML获取照片位姿和共视信息
├── 构建候选像对并计算质量评分
└── 建立邻接表用于三角形搜索

Phase 2: 三角形发现与选择
├── 对每对候选像对(i,j)搜索公共邻居k
├── 评估三角形(i,j,k)质量
├── 按质量排序所有候选三角形
└── 贪心选择三角形直到达到目标数量

Phase 3: 密集匹配对生成 (区段1)
├── 从选中的三角形提取所有边
├── 选择高质量子集作为密集匹配对
└── 确保每个照片有足够邻居

Phase 4: Refine对生成 (区段2)
├── 以密集匹配对为基础
├── 添加部分高共视但未被选中的边
└── 确保全局连通性和BA优化需求

Phase 5: 三段式输出
├── 区段1: 密集匹配对
├── 区段2: Refine对
└── 区段3: 三角形约束
```

### 2.3 三角形质量评估

```cpp
double evaluateTriangle(int i, int j, int k) {
    // 1. 三边共视度
    int covis_ij = getCovis(i, j);
    int covis_jk = getCovis(j, k);
    int covis_ik = getCovis(i, k);
    int minCovis = min({covis_ij, covis_jk, covis_ik});

    // 2. 弱边检查 - 任一边太弱则放弃
    if (minCovis < MIN_COVIS_THRESHOLD) return 0;

    // 3. 角度分布 (避免退化三角形)
    double angle_ij = getAngle(i, j);
    double angle_jk = getAngle(j, k);
    double angle_ik = getAngle(i, k);
    double meanAngle = (angle_ij + angle_jk + angle_ik) / 3.0;

    // 4. 综合评分
    double covisScore = (double)minCovis / MAX_COVIS;
    double angleScore = 1.0 - min(meanAngle / 45.0, 1.0);

    return 0.8 * covisScore + 0.2 * angleScore;
}
```

### 2.4 密集对选择策略

```cpp
// 从三角形边中选择密集匹配对
std::set<Edge> selectDensePairs(const std::vector<Triplet>& triplets) {
    // 1. 收集所有三角形边及其出现次数
    std::map<Edge, int> edgeCount;
    std::map<Edge, double> edgeScore;

    for (auto& tri : triplets) {
        for (auto& edge : tri.getEdges()) {
            edgeCount[edge]++;
            edgeScore[edge] = max(edgeScore[edge],
                                  getEdgeScore(edge.first, edge.second));
        }
    }

    // 2. 优先选择出现在多个三角形中的边
    std::vector<Edge> sortedEdges;
    for (auto& [edge, count] : edgeCount) {
        sortedEdges.push_back(edge);
    }

    // 按(出现次数, 评分)排序
    sort(sortedEdges.begin(), sortedEdges.end(),
         [&](const Edge& a, const Edge& b) {
             if (edgeCount[a] != edgeCount[b])
                 return edgeCount[a] > edgeCount[b];
             return edgeScore[a] > edgeScore[b];
         });

    // 3. 选择目标数量的边
    std::set<Edge> densePairs;
    int targetCount = (int)(triplets.size() * 1.8);  // 约为三角形数×1.8

    for (auto& edge : sortedEdges) {
        if (densePairs.size() >= targetCount) break;
        densePairs.insert(edge);
    }

    return densePairs;
}
```

### 2.5 目标参数

| 参数 | 目标值 | 说明 |
|------|--------|------|
| 三角形/照片比 | ~0.63-0.69 | 每个照片约参与0.6-0.7个三角形 |
| 密集对/三角形比 | ~1.8 | 密集对约为三角形数×1.8 |
| Refine/Dense比 | ~1.0-1.02 | Refine略多或相等 |
| 密集对覆盖率 | >95% | 密集对应被三角形覆盖 |

## 3. 数据结构设计

```cpp
// 三角形结构
struct Triplet {
    int id1, id2, id3;  // 顶点ID (已排序: id1 < id2 < id3)
    double quality;      // 三角形质量

    // 获取三条边
    std::vector<std::pair<int,int>> getEdges() const {
        return {{id1, id2}, {id1, id3}, {id2, id3}};
    }

    // 用于去重的比较
    bool operator<(const Triplet& other) const {
        if (id1 != other.id1) return id1 < other.id1;
        if (id2 != other.id2) return id2 < other.id2;
        return id3 < other.id3;
    }
};

// 选择结果
struct PairSelectionResult {
    std::vector<std::pair<int,int>> densePairs;   // 区段1: 密集匹配
    std::vector<std::pair<int,int>> refinePairs;  // 区段2: Refine
    std::vector<Triplet> triplets;                // 区段3: 三角形
};
```

## 4. 输出格式

```cpp
void writeThreeSectionPairs(const std::string& path,
                            const PairSelectionResult& result) {
    std::ofstream file(path);

    file << "#Graph of Views\n";
    file << "#Version 2\n";

    // 区段1: 密集匹配对
    file << result.densePairs.size() << "\n";
    for (const auto& p : result.densePairs) {
        file << p.first << " " << p.second << "\n";
    }

    // 区段2: Refine对
    file << result.refinePairs.size() << "\n";
    for (const auto& p : result.refinePairs) {
        file << p.first << " " << p.second << "\n";
    }

    // 区段3: 三角形
    file << result.triplets.size() << "\n";
    for (const auto& t : result.triplets) {
        file << t.id1 << " " << t.id2 << " " << t.id3 << "\n";
    }
}
```

## 5. 实现计划

### Phase 1: 三角形发现模块 (TripletFinder)
- [x] 构建共视邻接表
- [ ] 实现公共邻居搜索
- [ ] 实现三角形质量评估
- [ ] 贪心三角形选择

### Phase 2: 密集对提取模块 (DensePairExtractor)
- [ ] 从三角形提取边并统计出现次数
- [ ] 边排序和选择
- [ ] 连通性保证

### Phase 3: Refine对生成模块 (RefinePairGenerator)
- [ ] 以密集对为基础
- [ ] 添加高共视候选
- [ ] 确保BA优化需求

### Phase 4: 集成输出
- [ ] 三段式输出格式
- [ ] 与参考输出对比
- [ ] 参数调优

## 6. 预期效果

| 指标 | 当前 | 目标 |
|------|------|------|
| F1 (vs 全部) | ~73% | >80% |
| 密集对匹配率 | ~65% | >85% |
| 三角形匹配率 | N/A | >80% |
| 输出格式 | 单段 | 三段 |

## 7. 关键参数

```cpp
struct TriangleConfig {
    // 三角形选择
    double targetTripletRatio = 0.65;       // 目标三角形/照片比
    int minTriangleCovis = 30;              // 三角形最小边共视度
    double minTriangleQuality = 0.3;        // 最小三角形质量

    // 密集对选择
    double denseToTripletRatio = 1.8;       // 密集对/三角形比
    int minDenseCovis = 50;                 // 密集对最小共视度

    // Refine对选择
    double refineExpansion = 1.02;          // Refine扩展比例
    int minRefineCovis = 30;                // Refine最小共视度
};
```

---
*更新时间: 2025-01-07*
*版本: v2.0 - 三段式输出策略*
