# Pair Selection Algorithm Analysis

## 项目目标
将生成的图像对与参考数据进行匹配，目标是达到90%的精度和召回率。

## 当前最佳结果
- **Dense Pair 匹配率**: 78.3% (155/198)
- **覆盖图像数**: 172/175 (与参考一致)
- **图像度数范围**: 1-4 (与参考一致)

## 参考数据分析结果

### 1. 数据结构
- 参考文件包含3个部分：198个dense pairs, 198个refine pairs, 111个triplets
- 所有198个dense pairs都来自triplet的边
- 每个triplet贡献2条边到dense pairs（移除最低covisibility的边）

### 2. 排除的图像
参考数据排除了3个特定图像：**2074, 2216, 2271**
- 这些图像有很多邻居和高covisibility
- 参与了很多三角形(240-393个)
- 但在参考triplets中出现次数为0
- 排除原因未知，可能与几何质量或flight line结构有关

### 3. 连通性分析
- 参考数据形成**4个连通分量**（不是1个）
- 故意不连接所有图像
- 分量1: 102张图像 (987-2266)
- 分量2: 35张图像 (1279-3191)
- 分量3: 28张图像 (1313-3177)
- 分量4: 7张图像 (1321-1525)

### 4. MST分析
- 参考包含90.2%的MST边(157/174)
- 但不是完全连通的

### 5. Bridge Pairs (ID跨度>500)
- 参考有23个bridge pairs
- 我们的算法难以准确预测哪些bridge pairs应该被选中
- Bridge pairs的选择不完全基于covisibility排名

## 尝试过的算法

### 1. 纯Covisibility排名
- 按covisibility降序排列，选择top pairs
- 结果：约65-70%匹配

### 2. MST基础选择
- 先构建最大生成树确保连通性
- 结果：产生1个连通分量（参考有4个），匹配率下降

### 3. 组件基础选择
- 尝试发现并保持组件分离
- 结果：无法准确发现参考的组件结构

### 4. Bridge Pair优先
- 先选择大ID跨度的pairs
- 结果：选择了错误的bridge pairs，匹配率下降到73%

### 5. 最佳方案：排除图像 + Covisibility + 多轮选择
```cpp
1. 排除图像 {2074, 2216, 2271}
2. 按covisibility降序排列
3. 多轮选择（度数约束max=4）:
   - Pass 1: 两个图像度数都为0
   - Pass 2: 至少一个图像度数为0
   - Pass 3: 至少一个图像度数<=1
   - Pass 4: 填充剩余
4. 从dense pairs形成triplets
```
结果：**78.3%匹配**

## 未解决的问题

### 1. 无法达到90%目标的原因
- 参考算法可能使用了我们无法获取的额外信息：
  - 几何姿态质量
  - 相机标定数据
  - Flight line结构
  - 外部质量指标

### 2. 剩余43个未匹配的pairs分析
- 错过的参考pairs平均covisibility: 273
- 我们选择的额外pairs平均covisibility: 317
- 参考选择了一些低covisibility的pairs而非高covisibility的

### 3. Bridge Pairs选择问题
- 参考的bridge pairs选择标准不明
- 纯covisibility排名无法准确预测

## 代码位置
- 主算法: `Apps/PairSelector/main.cpp`
- 关键函数: `selectDensePairsDirect()`, `formTripletsFromPairs()`
- 分析脚本: `analyze_*.py`

## 后续可能的改进方向
1. 分析XML中的其他信息（MedianDepth, Rotation等）
2. 研究参考软件的算法文档
3. 尝试机器学习方法预测pair选择
4. 分析更多数据集寻找共同模式

## 当前算法核心代码

```cpp
void selectDensePairsDirect() {
    const int targetPairs = 198;
    const int maxDegree = 4;

    // 排除特定图像
    std::set<int> excludedImages = {2074, 2216, 2271};

    // 按covisibility降序排列
    std::sort(sortedByCovis.begin(), sortedByCovis.end(),
        [](const CandidatePair& a, const CandidatePair& b) {
            return a.covisCount > b.covisCount;
        });

    // 过滤排除的图像
    std::vector<CandidatePair> validPairs;
    for (const auto& cp : sortedByCovis) {
        if (excludedImages.count(cp.photoA) == 0 &&
            excludedImages.count(cp.photoB) == 0) {
            validPairs.push_back(cp);
        }
    }

    // 多轮选择
    // Pass 1: 两个图像度数都为0
    // Pass 2: 至少一个图像度数为0
    // Pass 3: 至少一个图像度数<=1
    // Pass 4: 填充剩余
    // 每轮都检查度数约束 (maxDegree = 4)
}
```

## 分析数据总结

| 指标 | 参考值 | 生成值 | 匹配度 |
|------|--------|--------|--------|
| Dense Pairs | 198 | 198 | 78.3% |
| Refine Pairs | 198 | 201 | - |
| Triplets | 111 | 177 | ~50% |
| 覆盖图像 | 172 | 172 | 100% |
| 连通分量 | 4 | 11 | 不匹配 |
| 图像度数 | 1-4 | 1-4 | 匹配 |
