# SOP：GBM / 肿瘤微环境单细胞分析工作流

> **目标**：让 AI（或分析者）拿到一份新的 10x 单细胞数据后，不用摸索参数/代码/配色，直接照本流程产出可分群、可注释、可做双靶点 CAR 评估的结果。
>
> 适用：GBM / 实体瘤微环境（含肿瘤 + 髓系 + T/NK）。流程通用，marker 表可按需替换。
>
> 配套：`code/`（可复用脚本模板）、`config/`（参数与配色集中管理）。

---

## 0. 一句话总流程

```
下载 → 合并 → QC → log1p/HVG/PCA → BBKNN 整合 → 第一轮聚类(低res)分大谱系
→ 免疫定向提取 → 谱系内独立细分(肿瘤/T/髓系/NK) → 合并最终注释
→ 验证(rank_genes+dotplot) → 靶点分析(pseudobulk/共表达/AND-gate安全性) → 画图
```

---

## 1. 环境与依赖

- Python 3.12+，venv
- 安装（网络不通时用清华镜像）：`scanpy anndata leidenalg bbknn h5py matplotlib pandas numpy scipy pybiomart`

```bash
pip install scanpy anndata leidenalg bbknn -i https://pypi.tuna.tsinghua.edu.cn/simple
```

---

## 2. 数据准备

### 2.1 下载
- GEO 数据集，优先下载 `*_RAW.tar`（部分数据集个体文件不能单独下，需整个 tar）
- 断点续传用 `.part` 方案（参考 `dl_gse278450.py`）
- 若访问 NCBI 超时 → 开 VPN

### 2.2 蛋白编码基因过滤（可选但推荐）
- 用 pybiomart 查 Ensembl GRCh38 `biotype=protein_coding` 的 symbol
- 保留蛋白编码基因，排除 lncRNA 等

---

## 3. 预处理（照做）

> 脚本：`code/sc_preprocess.py`。改 `__main__` 里的文件路径即可。

| 步骤 | 动作 | 默认参数 |
|------|------|---------|
| 合并 | 逐样本读入(基因×细胞)→转稀疏→vstack | obs 存 sample/batch |
| QC | 过滤 | min_genes=200, mito<10%, n_genes<5000, counts<40000, min_cells=3 |
| 标准化 | normalize_total → **log1p** → 存 raw | target_sum=1e4 |
| HVG | 高变基因 | 2000, seurat, batch_key=sample |
| PCA | 仅对 HVG scale 后 | scale_max=10, n_comps=50 |

⚠️ 内存：90k 细胞 × 万基因勿 `.toarray()`，全程 sparse；scale 只对 HVG。

---

## 4. 整合 + 第一轮聚类（分大谱系）

> 脚本：`code/cluster_annotate.py` 的 `cluster_round1`

| 步骤 | 动作 | 参数 |
|------|------|------|
| 整合 | **BBKNN** | batch_key=sample, n_pcs=30, neighbors_within_batch=10 |
| 降维 | UMAP（基于 BBKNN 图）| — |
| 聚类 | Leiden | **resolution=0.3**（低，只分大谱系）|

> 结果：9~15 个大 cluster。经验：肿瘤量大会主导，低分辨率避免把肿瘤内部切碎。

---

## 5. 第一轮谱系注释

| 大谱系 | 代表 marker |
|--------|------------|
| Tumor | EGFR, PDGFRA, OLIG2, GFAP, SOX2, VIM |
| Myeloid | CD14, CD68, LYZ, CSF1R, AIF1, C1QA |
| Oligo | MBP, PLP1, MOG, CLDN11 |
| T | CD3D, CD3E, CD3G, CD2 |
| NK | NKG7, GNLY, KLRD1 |
| Endo | VWF, CLDN5, PECAM1 |

做法：每个 cluster 平均 marker 表达打分 + **手动核对覆盖**（例如增殖 cluster 归入其母谱系）。

---

## 6. 免疫定向提取（肿瘤组织关键）

> GBM/肿瘤中，第一轮 cluster 往往只抓住少数免疫细胞，真实 T/NK 散落各处，需按 marker 阈值提取。

| 细胞 | 判据 |
|------|------|
| T | CD3D/CD3E/CD3G/CD2 任一 log>0.5，排除强髓系 |
| NK | NKG7/GNLY/KLRD1 任一 >0.5，且 CD3 阴，排除强髓系 |
| NKT | CD3+ 且强 NK+ |

⚠️ 坑：GBM 肿瘤细胞可能弱表达免疫 marker → T/NK 会混入肿瘤。若需干净免疫群，用 CNV（inferCNV/copyKAT）区分恶性/非恶性，或反向排除肿瘤标志（SOX2/EGFR/GFAP/AQP4/VIM）。

---

## 7. 谱系内独立细分（核心技巧）

> **不要在全数据上高分辨率聚类**（肿瘤主导会压扁免疫/髓系）。对每个谱系**单独**重算 HVG+PCA+UMAP+Leiden。

| 谱系 | 期望类别 | Leiden res |
|------|---------|-----------|
| Tumor | AC / MES / NPC / OPC（4 状态）| 0.4 |
| Myeloid | Mono/Mac / Microglia / Proliferating | 0.4 |
| T | CD8 / CD4 / Memory / Treg | 0.5 |
| NK | CD56bright / Cytotoxic / Activated / CD56dim | 0.5 |

> ⚠️ 髓系注释别用 HLA-DR 判 DC（所有巨噬都高表达），DC 用 FCER1A/CLEC9A/CD1C。

---

## 8. 合并最终注释 + 验证

### 合并
各谱系细分标签按 barcode 映射回全体细胞 → `final_celltype`（如 `Tumor-AC`、`Myeloid-Microglia`、`T-CD8`…）。

### 验证（必须，auto-marker 不可信）
- **rank_genes_groups**（无监督）：每群 top 差异基因是否与注释吻合
- **dotplot**：marker 特异性
- 检查可疑群：top 基因若含其它谱系 marker（如 T 群 top 基因是 SOX2/EGFR）→ 混入肿瘤

---

## 9. 靶点 / 双靶点 CAR 分析

> 脚本：`code/target_analysis.py`

### 9.1 各靶点表达在哪个群（用 pseudobulk CPM 消除 dropout）
- 原始 count 按亚群 pooled → CPM
- 判断：CD70→T？IL13RA2→肿瘤？PILRA→髓系？

### 9.2 AND-gate 安全性（双阳窗口落哪）
- 对靶点对 (g1,g2)：双阳细胞落在肿瘤 vs 非肿瘤的比例
- **窗口偏向肿瘤** → 安全偏好（如 IL13RA2+PILRA: 68% 肿瘤内）
- **窗口偏向非肿瘤**（如 CD70+PILRA 落髓系）→ AND 门会误伤

### 9.3 注意
- mRNA（单细胞）≠ 表面蛋白（CAR 真靶点）→ 需流式/蛋白验证
- 双阳率往往很低（AND 门"安全但低效"）——展示时如实说明

---

## 10. 画图规范（统一风格）

> 脚本：`code/plot_utils.py`，配色统一读 `config/palette.json`。

| 图 | 用哪个配色/色标 |
|----|----------------|
| 全局分群 UMAP | palette_A（谱系分组色）|
| 大谱系 | palette_B |
| 各谱系注释级分群 | palette_C（高区分度，子状态不同色）|
| 三基因表达 featureplot | **viridis + log1p**（0=不表达）|
| 跨谱系靶点热图 | **z-score**（跨基因可比）|
| marker dotplot | log1p |

### log1p vs z-score 速记
| | log1p | z-score |
|--|-------|---------|
| 问什么 | 基因在哪表达、多少 | 相对高低、谁更高 |
| 0 = | 不表达 | 平均 |
| 负值 | 无 | 有 |
| 跨基因比 | 否 | 是 |
| 用于 | featureplot/dotplot | 跨亚群热图/聚类 |

---

## 11. 常见坑（备忘）

1. **稀疏数据勿 toarray**（OOM）
2. **tight_layout 会忽略手动 colorbar** → 用 bbox_inches='tight'
3. **obs 数组只读** → 赋值先 `.copy()`
4. **GBM 肿瘤/胶质 marker 重叠** → auto-marker 需验证；必要时 CNV 分恶性
5. **单细胞低表达基因 dropout 严重** → 用 pseudobulk 还原丰度；表达图用 log1p
6. **z-score 画稀疏单细胞图会"看不到负值"**（0 值细胞聚 0 附近）→ 勿用
7. 换数据集：改 `config/params.json` + 各脚本路径 + marker 表
