"""
target_analysis.py - 靶点共表达/AND-gate 分析可复用模板
=========================================================
用于双靶点 CAR 设计评估。核心问题:
  1. 各靶基因(如 CD70/IL13RA2/PILRA) 在哪些细胞群表达? (pseudobulk CPM / log1p)
  2. 双靶点(AND-gate)是否在同一细胞共表达? 双阳窗口落在肿瘤还是正常? (安全性)

归一化准则:
  - log1p: 单细胞表达图 (0=不表达)
  - pseudobulk CPM: 消除 dropout 的群体真实丰度 (count pooled)
  - z-score: 跨谱系/跨基因相对比较 (热图)
"""
import os, sys, io
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')


# ============================================================
# 1. 各亚群 pseudobulk CPM (消除 dropout)
# ============================================================
def pseudobulk_cpm(adata, groupby, genes, use_raw_count=True,
                   min_cells=30):
    """
    adata: 原始 count 的 AnnData (X 为 count)
    groupby: obs 中分群列 (e.g. final_celltype / subtype)
    genes: 关注的基因列表
    返回 DataFrame: 行=亚群, 列=基因 CPM + n_cells
    """
    Xc = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    gpos = {g: i for i, g in enumerate(adata.var_names)}
    rows = []
    for grp in adata.obs[groupby].unique():
        mask = (adata.obs[groupby] == grp).values
        n = int(mask.sum())
        if n < min_cells:
            continue
        tot = np.asarray(Xc[mask].sum(axis=0)).ravel()
        lib = tot.sum()
        if lib == 0:
            continue
        cpm = tot / lib * 1e6
        rows.append({'subtype': grp, 'n_cells': n,
                     **{g: float(cpm[gpos[g]]) for g in genes if g in gpos}})
    return pd.DataFrame(rows).set_index('subtype')


# ============================================================
# 2. 单细胞共表达 (co-detection) + Fisher exact
# ============================================================
def coexpr_table(adata, g1, g2, groupby, th=0.5):
    """
    对每个亚群算 g1,g2 共表达: 双阳数、各自检出率、双阳率、Fisher p.
    基于 log 表达 (adata.raw) 二值化 >th.
    """
    raw = adata.raw.to_adata()
    Xc = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    gpos = {g: i for i, g in enumerate(raw.var_names)}

    def vec(g):
        return np.asarray(Xc.getcol(gpos[g]).toarray()).ravel() if g in gpos \
            else np.zeros(Xc.shape[0])

    xb, yb = vec(g1) > th, vec(g2) > th
    rows = []
    for grp in adata.obs[groupby].unique():
        m = (adata.obs[groupby] == grp).values
        n = int(m.sum())
        if n < 10: continue
        a, b, ab = int((xb & m).sum()), int((yb & m).sum()), int((xb & yb & m).sum())
        tbl = [[ab, a - ab], [b - ab, n - a - b + ab]]
        try:
            _, p = stats.fisher_exact(tbl)
        except Exception:
            p = np.nan
        jacc = ab / (a + b - ab) if (a + b - ab) > 0 else 0
        exp_ab = (a / n) * (b / n) * n
        enrich = ab / exp_ab if exp_ab > 0 else np.nan
        rows.append({'group': grp, 'n': n, 'only_g1': a - ab, 'only_g2': b - ab,
                     'both': ab, 'rate_g1': a / n, 'rate_g2': b / n,
                     'rate_both': ab / n, 'fisher_p': p, 'jaccard': jacc,
                     'enrich': enrich})
    return pd.DataFrame(rows)


# ============================================================
# 3. AND-gate 安全性: 双阳窗口落在肿瘤还是非肿瘤
# ============================================================
def andgate_safety(adata, g1, g2, tumor_mask, th=0.5):
    """
    tumor_mask: bool array, 1=肿瘤细胞
    返回: 双阳总数、肿瘤内双阳数/率、非肿瘤双阳数/率, 判断窗口偏向
    """
    raw = adata.raw.to_adata()
    Xc = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    gpos = {g: i for i, g in enumerate(raw.var_names)}

    def vec(g):
        return np.asarray(Xc.getcol(gpos[g]).toarray()).ravel()

    both = (vec(g1) > th) & (vec(g2) > th)
    non_tumor = ~tumor_mask
    n_tum, n_non = tumor_mask.sum(), non_tumor.sum()
    res = {
        'both_total': int(both.sum()),
        'both_tumor': int((both & tumor_mask).sum()),
        'both_nontumor': int((both & non_tumor).sum()),
        'tumor_rate': float((both & tumor_mask).sum()) / n_tum,
        'nontumor_rate': float((both & non_tumor).sum()) / n_non,
    }
    res['tumor_frac_of_both'] = res['both_tumor'] / max(res['both_total'], 1)
    res['window_bias'] = ('tumor' if res['tumor_rate'] > res['nontumor_rate']
                          else 'nontumor')
    return res


if __name__ == '__main__':
    print("导入用模块。实际运行参考 GSE278450 项目 andgate_safety / pseudobulk_subtype / coexpr_all")
