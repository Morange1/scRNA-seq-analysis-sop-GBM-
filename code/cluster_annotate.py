"""
cluster_annotate.py - 两轮聚类 + 注释 可复用模板
=================================================
流程: BBKNN 整合 -> 第一轮低分辨率分大谱系 -> 谱系内细分 -> marker 注释
适用于 GBM/肿瘤微环境单细胞(流程通用, marker 按 GBM 默认, 可替换)。

关键概念:
  - 两轮聚类: 第一轮低分辨率(0.3)分大谱系; 第二轮谱系内独立重算 HVG+PCA+Leiden
  - 谱系内独立重算很关键: 避免大群(肿瘤)主导 HVG, 压扁免疫/髓系内部结构
  - auto-marker 打分只做初步注释, 需配 rank_genes/dotplot 验证(见 verify)
"""
import os, sys, io, json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ---------------- 大谱系 / 亚群 marker 表 (GBM 默认, 可改) ----------------
LINEAGE_MARKERS = {
    'Tumor':   ['EGFR', 'PDGFRA', 'OLIG2', 'GFAP', 'SOX2', 'VIM'],
    'Myeloid': ['CD14', 'CD68', 'LYZ', 'CSF1R', 'AIF1', 'C1QA'],
    'Oligo':   ['MBP', 'PLP1', 'MOG', 'CLDN11', 'MAG'],
    'T':       ['CD3D', 'CD3E', 'CD3G', 'CD2'],
    'NK':      ['NKG7', 'GNLY', 'KLRD1', 'KLRF1'],
    'Endo':    ['VWF', 'CLDN5', 'PECAM1'],
}
TUMOR_STATE_MARKERS = {
    'NPC': ['SOX2', 'ASCL1', 'EGFR', 'DCX'],
    'OPC': ['PDGFRA', 'OLIG2', 'NKX2-2', 'GPR17'],
    'AC':  ['GFAP', 'GLUL', 'AQP4', 'SLC1A3', 'APOE'],
    'MES': ['CHI3L1', 'CD44', 'VIM', 'ANXA1', 'SERPINE1'],
}
MYE_MARKERS = {
    'Mono/Mac': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN', 'FCGR3A'],
    'Microglia': ['P2RY12', 'TMEM119', 'SALL1', 'GPR34', 'CST3'],
    'Proliferating': ['MKI67', 'TOP2A', 'CENPF'],
    'DC': ['FCER1A', 'CLEC9A', 'CD1C'],  # 注意: 高特异, 勿用 HLA-DR 误判
}
T_MARKERS = {
    'CD8': ['CD8A', 'CD8B', 'GZMK', 'CTSW'],
    'CD4': ['CD4', 'IL7R', 'LEF1', 'TCF7'],
    'Treg': ['FOXP3', 'CTLA4', 'IL2RA'],
    'Memory': ['SELL', 'CCR7', 'TCF7', 'IL7R'],
}
NK_MARKERS = {
    'CD56bright': ['NCAM1', 'SELL', 'IL7R', 'XCL1'],
    'CD56dim': ['FCGR3A', 'KLRD1', 'PRF1', 'GZMB'],
    'Cytotoxic': ['PRF1', 'GZMB', 'GNLY', 'NKG7'],
    'Activated': ['MKI67', 'CENPF', 'IFNG'],
}


# ---------------- 工具: 某子集内独立重聚类 ----------------
def recluster_subset(ad, subset_bool, res, n_pcs=30, hvg=2000, key='sub'):
    """对 ad 中 subset_bool=True 的细胞, 用各自 raw(log全基因) 独立重聚类."""
    sub = ad[subset_bool].copy()
    s = sub.raw.to_adata()
    s = s[sub.obs_names].copy()
    s.obs = sub.obs.copy()
    sc.pp.highly_variable_genes(s, n_top_genes=min(hvg, s.n_vars - 1))
    if not s.var['highly_variable'].any():
        s.var['highly_variable'] = True
    s = s[:, s.var['highly_variable']].copy()
    sc.pp.scale(s, max_value=10)
    ncomp = min(n_pcs, s.n_vars - 1)
    sc.tl.pca(s, n_comps=ncomp, svd_solver='arpack')
    sc.pp.neighbors(s, n_neighbors=15, n_pcs=ncomp)
    sc.tl.umap(s)
    sc.tl.leiden(s, resolution=res, flavor='igraph', n_iterations=2,
                 key_added=key)
    return s   # 含 key 亚群标签 + 局部 X_umap


# ---------------- 工具: marker 打分注释 ----------------
def annotate_by_marker(s, key, marker_map, label_out):
    """对 s 的每个亚群(key), 按 marker_map 打分取最高者注释 -> obs[label_out]."""
    raw = s.raw.to_adata() if s.raw is not None else s
    Xc = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    gpos = {g: i for i, g in enumerate(raw.var_names)}
    clusters = sorted(s.obs[key].unique(), key=int)
    ann = {}
    for c in clusters:
        mask = (s.obs[key] == c).values
        scores = {}
        for ct, genes in marker_map.items():
            genes = [g for g in genes if g in gpos]
            if not genes:
                scores[ct] = 0; continue
            scores[ct] = float(Xc[mask][:, [gpos[g] for g in genes]].mean())
        ann[c] = max(scores, key=scores.get)
    s.obs[label_out] = [ann[c] for c in s.obs[key]]
    return s, ann


# ---------------- 主流程 ----------------
def cluster_round1(adata, p):
    """
    第一轮: BBKNN 整合 + UMAP + Leiden(低分辨率) -> obs['leiden_global']
    需要 adata 已有 X_pca + raw. p: params 的 ['integration']+['leiden']
    """
    import bbknn
    bbknn.bbknn(adata, batch_key=p['integration']['batch_key'],
                n_pcs=p['integration']['n_pcs'],
                neighbors_within_batch=p['integration']['neighbors_within_batch'])
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=p['leiden']['round1_resolution'],
                 flavor='igraph', n_iterations=2, key_added='leiden_global')
    return adata


def annotate_lineage(adata):
    """第一轮 cluster 映射到大谱系(需 manual 修正表)."""
    raw = adata.raw.to_adata()
    Xc = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    gpos = {g: i for i, g in enumerate(raw.var_names)}
    clusters = sorted(adata.obs['leiden_global'].unique(), key=int)
    ann = {}
    for c in clusters:
        mask = (adata.obs['leiden_global'] == c).values
        scores = {}
        for lin, genes in LINEAGE_MARKERS.items():
            genes = [g for g in genes if g in gpos]
            if not genes: scores[lin] = 0; continue
            scores[lin] = float(Xc[mask][:, [gpos[g] for g in genes]].mean())
        ann[c] = max(scores, key=scores.get)
    adata.obs['lineage'] = [ann[c] for c in adata.obs['leiden_global']]
    return adata, ann


if __name__ == '__main__':
    print("导入用模块。实际运行参考 GSE278450 项目脚本 step2_bbknn_round1 + subcluster_*")
