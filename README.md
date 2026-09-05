# scRNA-seq GBM/TME 分析 skill 包

针对 **GBM / 实体瘤微环境单细胞分析**的可复用工作流规范与代码模板。
核心目的：下次让 AI（或分析者）处理类似 10x 单细胞数据时，**无需再摸索流程、参数、代码与配色**，直接照 SOP 产出可分群、可注释、可评估双靶点 CAR 靶点的结果。

> **【使用守则 · 给 AI】** 在开始分析前，请通读本 README 与 `config/params.json`，随后**先评估默认参数的"可用性"**（判断它们对当前数据集是否合理），**列出针对本次单细胞分析需要改动的参数及理由**，**告知使用者**；**经使用者批准后，方可修改参数并继续分析**。未经批准不得擅自调整参数。

## 目录结构

```
skill/
├── SOP.md              # 指令式操作流程（照做即可，含参数/坑/准则）
├── config/
│   ├── params.json     # 所有参数集中管理（QC/HVG/聚类/marker阈值/靶基因）
│   └── palette.json    # 统一配色方案（分群三套 + 表达 colormap）
├── code/
│   ├── sc_preprocess.py    # 预处理：合并/QC/normalize/log1p/HVG/PCA
│   ├── cluster_annotate.py # 两轮聚类 + marker 注释（BBKNN/谱系内细分）
│   ├── target_analysis.py  # 靶点分析：pseudobulk/共表达/AND-gate安全性
│   └── plot_utils.py       # 统一画图：分群/表达/热图/dotplot（读 palette.json）
├── examples/
│   └── GSE278450_GBM/      # 真实跑通的附例（参数 + 结果摘要 + 双靶点结论）
└── README.md
```

## 快速开始（新数据集）

1. 看 `SOP.md`（流程 + 参数 + 坑）
2. 改 `config/params.json`（数据路径、QC、resolution、靶基因）
3. 改各脚本 `__main__` 的输入文件路径
4. 依次运行：`sc_preprocess → cluster_annotate → target_analysis → plot_utils` 对应用例
5. 如需换细胞类型：替换各脚本顶部的 marker 表

## 核心设计要点

- **两轮聚类**：第一轮低分辨率分大谱系；第二轮**谱系内独立重算** HVG+PCA，避免大群(肿瘤)主导
- **配色统一**：分群(谱系分组/高区分度两套) + 表达(viridis+log1p)，集中存于 `palette.json`
- **log1p vs z-score 准则**：
  - 单细胞表达图用 **log1p**（0=不表达，直观）
  - 跨谱系/基因比较热图用 **z-score**（消除量级，可比）
- **AND-gate 双靶点评估**：pseudobulk CPM（消 dropout）+ Fisher 共表达 + 双阳窗口(肿瘤 vs 非肿瘤)安全性

## 参数与真实附例

### `config/params.json`（默认参数）
所有分析参数的**唯一权威来源**（QC 阈值 / normalize / HVG / 整合 / Leiden 分辨率 / marker 阈值 / 靶基因）。
换数据集时**优先改这里**，再改各脚本的数据路径。各脚本已按 `params.json` 结构编写，读同一份配置避免参数散落。

### `examples/GSE278450_GBM/`（真实跑通附例）
一个**完整跑通的参考案例**，展示本 skill 在真实 22 例 GBM 上"参数取值 → 结果产出 → 配色应用 → 双靶点结论"全过程：

| 文件 | 内容 |
|------|------|
| `params_actual.json` | GSE278450 真实参数取值 + 每阶段的真实细胞/基因数（91,336 → 90,609 → 17 类）|
| `results_summary.md` | 最终分群表、配色使用示例、log1p vs z-score 用法、双靶点 AND-gate 结论（IL13RA2+PILRA 更优）、复盘踩坑 |

> 用法建议：新数据集先用默认 `params.json` 跑，遇到问题对照附例的"参数取值与结果合理性"，看产出是否在合理范围。

## 参考数据集

本规范基于 **GSE278450**（Couturier et al. 2022, Neuro-Oncology, 22 例 GBM）分析沉淀。
原始处理脚本与结果在本地 `D:\GBM_neftel_project\01_scripts` 与 `02_results`（h5ad 数据不随仓库发布）。

## 环境依赖

```
python>=3.12
scanpy, anndata, leidenalg, bbknn, h5py, scipy, numpy, pandas, matplotlib, pybiomart
```
