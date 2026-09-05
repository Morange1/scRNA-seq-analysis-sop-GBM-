"""
sc_preprocess.py - 单细胞预处理可复用模板
===========================================
针对 10x 单细胞 RNA-seq（GBM/肿瘤微环境为主，流程通用）。

用法: 换数据集时, 只需:
  1. 修改 __main__ 里的 DATA 配置 (输入文件/输出目录)
  2. 运行, 得到 QC + normalize + log1p + HVG + PCA 后的 AnnData
  3. 进入下一步 cluster_pipeline.py

依赖: scanpy, anndata, scipy, pandas, numpy
归一化准则: 本模块做 log1p (表达量绝对视图, 0=不表达)。
          z-score 仅在聚类输入/跨谱系热图时用(见 annotate/plot 模块)。
"""
import os
import sys
import io
import json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')


# ============================================================
# 1. 从多个样本的 10x matrix (基因x细胞, 每样本一个文件) 合并
# ============================================================
def merge_10x_samples(raw_files, sample_names, keep_coding=True,
                      coding_symbols=None):
    """
    raw_files: list of 10x 输出文件路径 (tsv.gz, 基因x细胞) 或 mtx 三元组目录
    sample_names: 对应样本名 (写入 obs['sample'], batch)
    返回合并后 AnnData (原始 count, CSR)
    """
    all_X, all_obs = [], []
    gene_list = None
    for i, fp in enumerate(raw_files):
        print(f'  读入 {sample_names[i]} ...')
        df = pd.read_csv(fp, sep='\t', index_col=0, engine='c')
        genes = df.index.astype(str).tolist()
        if gene_list is None:
            gene_list = genes
        else:
            assert genes == gene_list, '各样本基因列表不一致!'
        X = sp.csr_matrix(df.values.T)  # 细胞 x 基因
        cells = df.columns.astype(str).tolist()
        all_X.append(X)
        all_obs.append(pd.DataFrame({'sample': sample_names[i],
                                     'batch': sample_names[i]}, index=cells))
    Xbig = sp.vstack(all_X).tocsr()
    obs = pd.concat(all_obs)
    adata = sc.AnnData(X=Xbig, obs=obs)
    adata.var_names = gene_list
    # 可选: 过滤到蛋白编码基因
    if keep_coding and coding_symbols is not None:
        coding = set(coding_symbols)
        keep = [g for g in adata.var_names if g in coding]
        adata = adata[:, keep].copy()
    return adata


# ============================================================
# 2. QC (按 params.json 中的阈值)
# ============================================================
def qc_filter(adata, p):
    """p: params dict 的 ['qc'] 部分. 返回过滤后的 adata + 统计"""
    mito_prefix = p.get('mito_prefix', 'MT-')
    mito = adata.var_names.str.startswith(mito_prefix)
    if mito.sum() > 0:
        tot = np.asarray(adata.X.sum(axis=1)).ravel()
        tot[tot == 0] = 1
        adata.obs['pct_mito'] = \
            np.asarray(adata[:, mito].X.sum(axis=1)).ravel() / tot * 100
    else:
        adata.obs['pct_mito'] = 0.0
    adata.obs['n_genes'] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata.obs['total_counts'] = np.asarray(adata.X.sum(axis=1)).ravel()
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=p['min_genes'])
    adata = adata[adata.obs['pct_mito'] < p['max_pct_mito'], :].copy()
    adata = adata[adata.obs['n_genes'] < p['max_genes'], :].copy()
    adata = adata[adata.obs['total_counts'] < p['max_total_counts'], :].copy()
    sc.pp.filter_genes(adata, min_cells=p['min_cells'])
    print(f'  QC: {n0} -> {adata.n_obs} 细胞')
    return adata


# ============================================================
# 3. 标准化 + HVG + PCA (全流程一步)
# ============================================================
def normalize_hvg_pca(adata, p):
    """
    标准单细胞处理: normalize_total(log1p) -> raw -> HVG -> scale -> PCA
    p: params 的 ['normalize']+['hvg']+['pca']
    返回: adata (X 保持 log1p 全基因, raw=log1p 全基因, obsm['X_pca'])
    """
    # log1p (绝对表达视图; z-score 仅在 cluster/热图时用)
    sc.pp.normalize_total(adata, target_sum=p['normalize']['target_sum'])
    sc.pp.log1p(adata)
    adata.raw = adata.copy()          # log1p 全基因, 供 marker 注释/表达图
    # HVG
    sc.pp.highly_variable_genes(
        adata, n_top_genes=p['hvg']['n_top_genes'],
        flavor=p['hvg']['flavor'], batch_key=p['hvg']['batch_key'])
    # scale 只对 HVG 做 (避免全基因密集化 OOM)
    adata = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata, max_value=p['pca']['scale_max'])
    # PCA
    ncomp = min(p['pca']['n_comps'], adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=ncomp, svd_solver='arpack')
    return adata


# ============================================================
# 4. 免疫/细胞类型 marker 定向提取 (可选)
# ============================================================
def extract_immune(adata, t_markers=None, nk_markers=None,
                   mye_exclude=None, th=0.5):
    """
    从全体细胞中按 marker 阈值定向提取免疫细胞(T/NK)。
    注意: 肿瘤组织(尤其GBM)中肿瘤细胞可能弱表达免疫marker, 结果需复核。
    """
    raw = adata.raw.to_adata()
    Xc = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    gpos = {g: i for i, g in enumerate(raw.var_names)}

    def colpos(g):
        if g not in gpos:
            return np.zeros(Xc.shape[0], bool)
        v = np.asarray(Xc.getcol(gpos[g]).toarray()).ravel()
        return v > th

    t_markers = t_markers or ['CD3D', 'CD3E', 'CD3G', 'CD2']
    nk_markers = nk_markers or ['NKG7', 'GNLY', 'KLRD1']
    mye_exclude = mye_exclude or ['CD14', 'LYZ', 'CSF1R', 'AIF1', 'P2RY12']

    cd3_pos = np.logical_or.reduce([colpos(g) for g in t_markers])
    nk_pos = np.logical_or.reduce([colpos(g) for g in nk_markers])
    mye_hi = np.logical_or.reduce([colpos(g, th=2.0) for g in mye_exclude])

    t_mask = cd3_pos & (~mye_hi)
    nk_mask = nk_pos & (~cd3_pos) & (~mye_hi)   # NK: 非T
    nkt_mask = cd3_pos & (np.logical_or.reduce([colpos(g) for g in ['NKG7','GNLY']]))
    nkt_mask = nkt_mask & (~t_mask)  # 简化示意, 可按需调整
    adata.obs['t_mask'] = t_mask
    adata.obs['nk_mask'] = nk_mask
    return adata


if __name__ == '__main__':
    # ===== 按需修改: 数据路径 =====
    DATA = {
        'raw_files': [],
        'sample_names': [],
        'coding_file': None,
        'params_json': r'..\config\params.json',
        'out_h5ad': 'preprocessed.h5ad',
    }
    with open(DATA['params_json'], encoding='utf-8') as f:
        P = json.load(f)
    adata = merge_10x_samples(DATA['raw_files'], DATA['sample_names'])
    adata = qc_filter(adata, P['qc'])
    adata = normalize_hvg_pca(adata, P)
    adata.write(DATA['out_h5ad'])
    print('预处理完成 ->', DATA['out_h5ad'])
