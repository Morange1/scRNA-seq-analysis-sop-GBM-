# 附例：GSE278450 GBM 分析结果摘要

> 这是本 skill 在 **GSE278450**（Couturier et al. 2022, Neuro-Oncology；22 例 GBM）上真实跑通的参考案例。展示参数产出、配色应用、双靶点结论。新数据集可对照此例预估结果与检查合理性。

## 1. 数据概览（真实）

| 阶段 | 数值 |
|------|------|
| 输入 | 22 样本, 91,336 细胞 × 36,116 基因 |
| 蛋白编码过滤后 | 91,336 × 18,068 |
| QC 后 | **90,609 × 16,684**（727 cells mito>10% 被滤）|
| 免疫定向提取 | T=557, NK=1002, NKT=429 |
| 最终 final_celltype | 17 类, 90,609 细胞 |

## 2. 参数（见 params_actual.json）
关键：BBKNN(batch=sample, n_pcs=30)；round1 Leiden res=0.3；谱系内 T/NK res=0.5、肿瘤/髓系 res=0.4。

## 3. 最终细胞分群（真实产出）

| 大谱系 | 亚群 | 细胞数 |
|--------|------|--------|
| Tumor | AC / MES / NPC / OPC | 20,704 / 14,056 / 7,975 / 2,664 |
| Myeloid | Mono/Mac / Microglia / Proliferating | 21,892 / 10,885 / 1,580 |
| T | CD8 / CD4 / Memory / Treg | 317 / 163 / 64 / 13 |
| NK | CD56bright / Cytotoxic / Activated / dim | 508 / 363 / 73 / 58 |
| 其他 | NKT / Oligo | 429 / 8,865 |

## 4. 配色应用示例（读取 config/palette.json）

- **全局 17 类 UMAP** → palette_A（Tumor 红橙 / Myeloid 蓝 / T 绿 / NK 紫 / NKT 灰 / Oligo 黄）
- **大谱系 UMAP** → palette_B（Tumor 红 `#e41a1c` / Myeloid 蓝 / T 绿 / NK 紫 / Oligo 橙）
- **各谱系注释级分群** → palette_C（高区分度：如肿瘤 AC红/MES橙/NPC蓝/OPC绿）
- **三靶点表达 featureplot** → viridis + log1p（0=不表达深紫，表达渐绿→黄）

## 5. log1p vs z-score 在本例的体现

| 用途 | 本例做法 |
|------|---------|
| CD70/IL13RA2/PILRA 单细胞 UMAP 表达 | **log1p**（直观看到 IL13RA2 在肿瘤、PILRA 在髓系、CD70 零星）|
| 跨 5 大谱系比较三靶点相对高低 | **z-score**（三基因量级差大，z-score 后一张热图可比）|
| 经验 | z-score 画稀疏单细胞表达图会出现"大量细胞近 0、看不到负值"→ 改用 log1p |

## 6. 双靶点 CAR 结论（本例核心产出）

### 三靶点分布（pseudobulk CPM 佐证）
- **CD70** → 主表达于 T 细胞（T-CD8 CPM 39.5），肿瘤低
- **IL13RA2** → 主表达于肿瘤（MES CPM 93, NPC 47），肿瘤特异
- **PILRA** → 主表达于髓系（Microglia CPM 135, Mono/Mac 98）

### AND-gate 安全性（双阳窗口落肿瘤 or 非肿瘤）
| 配对 | 双阳窗口 | 肿瘤双阳率 | 非肿瘤双阳率 | 评价 |
|------|---------|-----------|-------------|------|
| IL13RA2 + PILRA | 68.3% 在肿瘤 | 0.300% | 0.137% | ✅ 安全窗口偏肿瘤 |
| CD70 + PILRA | 5.6% 在肿瘤 | 0.029% | 0.484% | ❌ 窗口偏髓系(误伤) |

**结论**：IL13RA2 + PILRA 是更优的双靶点 AND-gate 组合（覆盖肿瘤+髓系免疫抑制，窗口偏好肿瘤）。注意双阳率整体低（AND 门"安全但低效"），且 mRNA 需蛋白/流式层面验证。

## 7. 复盘与坑（本例实际遇到的）

1. **肿瘤细胞量大主导** → 必须谱系内独立重聚类，否则免疫被压扁
2. **髓系 DC 误判**：HLA-DR 全场巨噬高表达，误把巨噬标成 DC → 改用 FCER1A/CLEC9A/CD1C
3. **T/NK 混入肿瘤**：rank_genes 验证发现 T-CD4/Memory、NK-CD56bright/dim top 基因含 SOX2/EGFR 等肿瘤标志（GBM 冷肿瘤 + marker 重叠）。决定保留混杂（不影响双靶点核心结论）
4. **稀疏数据内存**：勿 `.toarray()` 全矩阵；只对 HVG scale
5. **obs 只读**：scanpy obs 底层数组只读，赋值需 `.copy()`

## 8. 产出文件参考（本地，不入库）
- 数据 h5ad：`00_data/GSE278450/gse278450_*.h5ad`（merged/pca/round1/final/各亚群）
- 图：`02_results/figures_278450/`（分群/表达/热图/dotplot）
- 表：`02_results/tables_278450/`（pseudobulk/coexpr/detect rate 等）
