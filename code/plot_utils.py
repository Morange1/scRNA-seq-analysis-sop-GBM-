"""
plot_utils.py - 统一配色 + 标准图函数
======================================
所有单细胞图的风格在这里统一:
  - 分群: palette_A (全局谱系) / palette_C (各谱系高区分度)
  - 表达: viridis, log1p (0=不表达)
  - 布局: 固定点大小/图例/colorbar 位置, 避免 tight_layout 与手动 colorbar 冲突

关键经验:
  - 稀疏数据(大量0)画单细胞表达用 log1p (直观: 0=不表达)
  - z-score 只用于跨亚群/基因热图 (量级不同需标准化)
  - tight_layout() 会重排子图但忽略手动 add_axes 的 colorbar -> 用 bbox_inches='tight'
"""
import os, sys, io, json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

PALETTE_PATH = os.path.join(os.path.dirname(__file__), '..', 'config', 'palette.json')
with open(PALETTE_PATH, encoding='utf-8') as f:
    PAL = json.load(f)

PAL_A = PAL['palette_A_lineage']       # 全局 17 类
PAL_B = PAL['palette_B_major']         # 大谱系
PAL_C = PAL['palette_C_subtype_highcontrast']  # 各谱系高区分度
CMAP_EXPR = PAL['expression_cmap']['cmap']    # 'viridis'


def _save(fig, out):
    # 不用 tight_layout, 避免与手动 colorbar 冲突
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print('  保存:', out)


# ============================================================
# 1. 全局分群 UMAP (谱系配色A)
# ============================================================
def plot_global_umap(adata, color='final_celltype', out='fig_subcluster_umap.png'):
    umap = adata.obsm['X_umap']
    if color == 'final_celltype':
        order = adata.obs[color].value_counts().index.tolist()
        # 每类颜色, 缺失用灰
        def col(ct): return PAL_A.get(ct, '#bbbbbb')
        fig, ax = plt.subplots(figsize=(10, 9))
        for ct in order:
            m = (adata.obs[color] == ct).values
            ax.scatter(umap[m, 0], umap[m, 1], color=col(ct), s=3, alpha=0.9,
                       rasterized=True)
        ax.legend(handles=[Patch(facecolor=col(ct), label=ct)
                           for ct in order], fontsize=7,
                  loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)
    else:
        sc.pl.umap(adata, color=color, show=False, ax=plt.gca(),
                   legend_loc='right margin')
    ax = plt.gca()
    ax.set_title(f'{color} ({adata.n_obs})', fontsize=13)
    ax.set_xticks([]); ax.set_yticks([])
    _save(plt.gcf(), out)


# ============================================================
# 2. 三基因 log1p 表达图 (全局或局部 UMAP)
# ============================================================
def plot_target_expression(adata, genes, umap_key='X_umap',
                           title='Target expression (log1p)',
                           out='fig_expr.png'):
    """
    adata.raw 需为 log1p 全基因. 每基因一图, viridis, 0=不表达(深紫).
    """
    raw = adata.raw.to_adata()
    Xc = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    gpos = {g: i for i, g in enumerate(raw.var_names)}
    umap = adata.obsm[umap_key]
    fig, axes = plt.subplots(1, len(genes), figsize=(6 * len(genes), 6))
    if len(genes) == 1:
        axes = [axes]
    for ax, g in zip(axes, genes):
        if g not in gpos:
            ax.set_title(f'{g} (NA)'); continue
        e = np.asarray(Xc.getcol(gpos[g]).toarray()).ravel()
        order = np.argsort(e)
        sc_ = ax.scatter(umap[order, 0], umap[order, 1], c=e[order],
                         cmap=CMAP_EXPR, s=2, alpha=0.9, rasterized=True)
        ax.set_title(f'{g} (log1p)', fontsize=14)
        ax.set_xticks([]); ax.set_yticks([])
    # colorbar 右侧, 不挤占
    pos = axes[-1].get_position()
    cax = fig.add_axes([pos.x1 + 0.015, pos.y0, 0.012, pos.height])
    cb = fig.colorbar(sc_, cax=cax); cb.set_label('log1p expression')
    fig.suptitle(title, fontsize=14)
    _save(fig, out)


# ============================================================
# 3. 单谱系高区分度分群 (palette C)
# ============================================================
def plot_lineage_cat(adata, label, palette, out='fig_cat.png'):
    umap = adata.obsm['X_umap']
    fig, ax = plt.subplots(figsize=(8, 7))
    for cat, c in palette.items():
        m = (adata.obs[label] == cat).values
        ax.scatter(umap[m, 0], umap[m, 1], color=c, s=3, alpha=0.9, label=cat,
                   rasterized=True)
    ax.set_title(f'{label} (n={adata.n_obs})', fontsize=13)
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(fontsize=9, loc='upper right')
    _save(fig, out)


# ============================================================
# 4. z-score 跨谱系热图 (基因级 z-score)
# ============================================================
def zscore_heatmap(cpm_df, genes, out='fig_zscore.png'):
    """
    cpm_df: 行=亚群, 列=基因 (CPM). 对每基因跨亚群 z-score.
    用于跨谱系比较 (不同基因量级需标准化). 单细胞表达图不要用z-score.
    """
    Z = (cpm_df[genes] - cpm_df[genes].mean(0)) / (cpm_df[genes].std(0) + 1e-9)
    fig, ax = plt.subplots(figsize=(max(3, len(genes) * 1.8),
                                    max(3, Z.shape[0] * 0.7)))
    im = ax.imshow(Z.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes)
    ax.set_yticks(range(Z.shape[0])); ax.set_yticklabels(Z.index)
    for i in range(Z.shape[0]):
        for j in range(len(genes)):
            ax.text(j, i, f'{Z.values[i, j]:.1f}', ha='center', va='center',
                    fontsize=8, color='white' if abs(Z.values[i, j]) > 1.2 else 'black')
    plt.colorbar(im, ax=ax, shrink=0.8, label='z-score')
    ax.set_title('z-score by lineage (pseudobulk CPM)', fontsize=12)
    _save(fig, out)


# ============================================================
# 5. dotplot marker 特异性 (验证注释)
# ============================================================
def marker_dotplot(adata, markers, groupby, out='fig_dotplot.png'):
    mk = [g for g in markers if g in (adata.raw.var_names if adata.raw else adata.var_names)]
    sc.pl.dotplot(adata, mk, groupby=groupby, show=False, standard_scale='var')
    _save(plt.gcf(), out)
