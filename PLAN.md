# 像对选择算法实现计划

## 项目目标

基于逆向工程分析结果，用C++实现像对选择算法，输入空三成果XML，输出pairs.txt。

---

## 一、核心算法设计

### 1.1 算法流程图

```
┌─────────────────┐
│  解析XML文件    │
│  - 照片位姿     │
│  - 连接点共视   │
└────────┬────────┘
         ▼
┌─────────────────┐
│  构建候选像对   │
│  - 共视点数>0   │
│  - 计算四因素   │
└────────┬────────┘
         ▼
┌─────────────────┐
│  综合评分       │
│  Score = f(O,A,G,D)│
└────────┬────────┘
         ▼
┌─────────────────┐
│  智能选择       │
│  - Top-K初选    │
│  - 冗余检测     │
│  - 三角形保证   │
└────────┬────────┘
         ▼
┌─────────────────┐
│  输出pairs.txt  │
└─────────────────┘
```

### 1.2 评分函数设计

基于分析结果，评分函数如下：

```cpp
double calculateScore(const PairFeatures& pair) {
    // 权重分配 (基于因素重要性排序)
    const double w_overlap = 0.50;    // 重叠度 - 最重要
    const double w_angle = 0.30;      // 交会角 - 非常重要
    const double w_gsd = 0.10;        // GSD比 - 较重要
    const double w_dir = 0.10;        // 方向相似度 - 较重要

    // 归一化处理
    double overlap_score = normalize_overlap(pair.covisCount);
    double angle_score = 1.0 - std::min(pair.convergenceAngle / 30.0, 1.0);
    double gsd_score = 1.0 - std::min(std::abs(pair.gsdRatio - 1.0) / 0.5, 1.0);
    double dir_score = pair.directionSimilarity;

    return w_overlap * overlap_score +
           w_angle * angle_score +
           w_gsd * gsd_score +
           w_dir * dir_score;
}
```

### 1.3 创新点实现

| 创新点 | 实现方法 |
|--------|---------|
| 智能冗余避免 | BFS检测2-hop连通性 |
| 非严格Top-K | 综合评分后按分数排序选择 |
| 三角形闭合 | 检查共同邻居数，优先保留能形成三角形的边 |
| 精确冗余控制 | 目标边数 ≈ 2 × (节点数-1) |
| 二级排序 | 同分数时优先选择小交会角 |

---

## 二、数据结构设计

### 2.1 核心数据结构

```cpp
// 照片信息
struct Photo {
    int id;
    double center[3];           // 相机中心坐标
    double rotation[9];         // 旋转矩阵 (M_00..M_22)
    double medianDepth;         // 中值深度
    int photogroup;             // 所属照片组
};

// 候选像对
struct CandidatePair {
    int photoA, photoB;
    int covisCount;             // 共视点数 (重叠度)
    double convergenceAngle;    // 交会角
    double gsdRatio;            // GSD比
    double directionSimilarity; // 方向相似度
    double score;               // 综合评分
};

// 图结构 (用于冗余检测)
class PairGraph {
    std::unordered_map<int, std::vector<int>> adjacency;

    bool is2HopConnected(int a, int b);
    int countCommonNeighbors(int a, int b);
    void addEdge(int a, int b);
};
```

### 2.2 XML解析模块

```cpp
class XMLParser {
public:
    bool parse(const std::string& xmlPath);

    const std::map<int, Photo>& getPhotos() const;
    const std::map<std::pair<int,int>, int>& getCovisibility() const;

private:
    void parsePhoto(const std::string& block);
    void parseTiePoint(const std::string& block);
};
```

---

## 三、实现步骤

### Step 1: XML解析器 (xml_parser.h/cpp)
- [ ] 解析Photo节点提取位姿信息
- [ ] 解析TiePoint节点统计共视关系
- [ ] 处理大文件（流式解析，避免全量加载）

### Step 2: 特征计算模块 (pair_features.h/cpp)
- [ ] 计算光轴方向 (从旋转矩阵)
- [ ] 计算交会角 (两光轴夹角)
- [ ] 计算GSD比 (通过中值深度)
- [ ] 计算方向相似度 (光轴余弦)

### Step 3: 评分与选择模块 (pair_selector.h/cpp)
- [ ] 综合评分函数
- [ ] Top-K初选
- [ ] 冗余检测 (BFS 2-hop)
- [ ] 三角形保证逻辑
- [ ] 最终选择

### Step 4: 输出模块 (output_writer.h/cpp)
- [ ] 生成pairs.txt格式
- [ ] 支持二列和三列格式

### Step 5: 测试验证
- [ ] 对4个数据集运行
- [ ] 与原始pairs.txt对比
- [ ] 统计匹配率和差异分析

---

## 四、项目结构

```
PairSelection/
├── PLAN.md                 # 本计划文件
├── CMakeLists.txt          # CMake构建配置
├── src/
│   ├── main.cpp            # 主程序入口
│   ├── xml_parser.h/cpp    # XML解析
│   ├── pair_features.h/cpp # 特征计算
│   ├── pair_selector.h/cpp # 选择算法
│   ├── pair_graph.h/cpp    # 图结构 (冗余检测)
│   └── output_writer.h/cpp # 输出模块
├── Datas/                  # 测试数据
│   ├── 1/
│   ├── 2/
│   ├── 3/
│   └── 4/
└── output/                 # 输出结果
```

---

## 五、参数配置

```cpp
struct Config {
    // 评分权重
    double weight_overlap = 0.50;
    double weight_angle = 0.30;
    double weight_gsd = 0.10;
    double weight_direction = 0.10;

    // 选择参数
    int minNeighborsPerPhoto = 2;   // 每张照片最少邻居数
    int maxNeighborsPerPhoto = 8;   // 每张照片最多邻居数
    double targetRedundancyRatio = 2.0;  // 目标冗余比 (相对于MST)

    // 质量阈值
    double maxConvergenceAngle = 30.0;  // 最大交会角 (度)
    double maxGsdRatio = 1.5;           // 最大GSD比
    int minCovisCount = 10;             // 最小共视点数
};
```

---

## 六、验证指标

### 6.1 定量指标

| 指标 | 目标值 | 验证方法 |
|------|--------|---------|
| 选择率 | ~3-6% | 选中数/候选数 |
| 平均度 | ~4 | 总边数×2/节点数 |
| 冗余比 | 1.8-2.2x | 边数/MST边数 |
| 聚类系数 | >0.4 | 计算平均聚类系数 |
| 连通性 | 100% | 验证图是否连通 |

### 6.2 与原始结果对比

```cpp
struct ComparisonResult {
    int originalPairs;      // 原始选中数
    int generatedPairs;     // 生成选中数
    int matchingPairs;      // 完全匹配数
    double precision;       // 精确率
    double recall;          // 召回率
    double f1Score;         // F1分数
};
```

---

## 七、时间估算

| 阶段 | 任务 | 文件 |
|------|------|------|
| 1 | XML解析器实现 | xml_parser.h/cpp |
| 2 | 特征计算实现 | pair_features.h/cpp |
| 3 | 图结构与冗余检测 | pair_graph.h/cpp |
| 4 | 选择算法实现 | pair_selector.h/cpp |
| 5 | 主程序与输出 | main.cpp, output_writer.cpp |
| 6 | 测试与验证 | 对比分析 |

---

## 八、风险与对策

| 风险 | 对策 |
|------|------|
| XML解析大文件内存不足 | 流式解析，逐块处理 |
| 冗余检测BFS性能问题 | 使用邻接表+visited数组优化 |
| 评分权重不准确 | 提供配置文件可调参 |
| 三角形保证过于严格 | 设置最小共同邻居阈值 |

---

## 九、下一步行动

1. **立即开始**: 创建项目结构和CMakeLists.txt
2. **优先实现**: XML解析器（核心依赖）
3. **迭代验证**: 每完成一个模块即运行测试
4. **最终输出**: 生成pairs.txt并对比验证

---

*计划创建时间: 2026-01-06*
*基于4个数据集的逆向工程分析*
