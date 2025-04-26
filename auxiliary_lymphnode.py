import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.stats import t
from statsmodels.stats.multitest import multipletests


def infer_celltype_activity(adata):
    A = adata.obsm['celltype_major'].to_numpy()
    b = adata.to_df().to_numpy()
    cov = np.dot(b.T - b.mean(), A - A.mean(axis=0)) / (b.shape[0]-1)
    ssd = np.std(A, axis=0, ddof=1) * np.std(b, axis=0, ddof=1).reshape(-1, 1)
    r = cov / ssd

    n_samples = b.shape[1]
    n_features, n_fsets = A.shape
    df = n_features - 2
    es = r * np.sqrt(df / ((1.0 - r + 1.0e-16)*(1.0 + r + 1.0e-16)))

    pv = t.sf(abs(es), df) * 2

    cts = adata.obsm['celltype_major'].columns
    tfs = adata.to_df().columns

    estimate = pd.DataFrame(es, index=tfs, columns=cts)
    estimate.name = 'coef'
    pvals = pd.DataFrame(pv, index=tfs, columns=cts)
    pvals.name = 'pvals'
    
    df_1 = estimate.stack().reset_index()
    df_1.columns = ['tf', 'ct', 'coef']
    
    df_2 = pvals.stack().reset_index()
    df_2.columns = ['tf', 'ct', 'pval']
    
    df_ct_tf = pd.merge(df_1, df_2, on=['tf', 'ct'])
    df_ct_tf["p_adj"] = multipletests(df_ct_tf['pval'], alpha=0.01, method="fdr_bh")[1]
    df_ct_tf["neg_log_p_adj"] = -np.log10(df_ct_tf["p_adj"]+1e-100)
    return df_ct_tf


def find_edges(adata, r=10, alpha=20, only_outer=True):
    points = np.asarray(adata[adata.obs['germinal_center']=='GC'].obsm['spatial'] * adata.uns['spatial']["V1_Human_Lymph_Node"]['scalefactors']['tissue_hires_scalef'])
    points = np.vstack((points+[-r,r], points+[-r,-r], points+[r,r], points+[r,-r]))
    assert points.shape[0] > 3, "Need at least four points"
    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))
    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return points, edges


def merge_celltypes(adata):
    df_celltype = pd.DataFrame(0, index=adata.obsm['celltype'].index, 
                          columns = ['B_Cycling', 'B_GC', 'B_IFN', 'B_activated', 'B_mem', 'B_naive', 'B_plasma', 'B_preGC', 
                                     'DC', 'Endo', 'FDC', 'ILC', 'Macrophages', 'Mast', 'Monocytes', 'NK', 'NKT', 
                                     'T_CD4+', 'T_CD8+', 'T_Treg', 'T_TIM3+', 'T_TfR', 'VSMC'])
    for ct in ['B_Cycling', 'B_IFN', 'B_activated', 'B_mem', 'B_naive', 'B_plasma', 'B_preGC', 
               'Endo', 'FDC', 'ILC', 'Mast', 'Monocytes', 'NK', 'NKT', 'T_Treg', 'T_TIM3+', 'T_TfR', 'VSMC']:
        df_celltype[ct] = adata.obsm['celltype'][ct]

    for ct in ['B_GC_DZ', 'B_GC_LZ', 'B_GC_prePB']:
        df_celltype['B_GC'] += adata.obsm['celltype'][ct]

    for ct in ['DC_CCR7+', 'DC_cDC1', 'DC_cDC2', 'DC_pDC']:
        df_celltype['DC'] += adata.obsm['celltype'][ct]

    for ct in ['Macrophages_M1', 'Macrophages_M2']:
        df_celltype['Macrophages'] += adata.obsm['celltype'][ct]

    for ct in ['T_CD4+', 'T_CD4+_TfH', 'T_CD4+_TfH_GC', 'T_CD4+_naive']:
        df_celltype['T_CD4+'] += adata.obsm['celltype'][ct]

    for ct in ['T_CD8+_CD161+', 'T_CD8+_cytotoxic', 'T_CD8+_naive']:
        df_celltype['T_CD8+'] += adata.obsm['celltype'][ct]
    return df_celltype


# def make_ct_tf_dataframe(adata, celltype_label='celltype'):
#     df_ct_tf = []
#     A = adata.obsm[celltype_label].copy()
#     Y = adata.to_df()
#     Y = Y-Y.mean()
#     Y = Y/Y.std()
#     for tf in adata.var_names:
#         y = Y[tf]
#         results = sm.OLS(y,A).fit()
#         for ct in A.columns:
#             p = results.pvalues[ct]
#             coef = results.params[ct]
#             se = results.bse[ct]
#             df_ct_tf.append([tf, ct, coef, p, se, results.rsquared])

#     df_ct_tf = pd.DataFrame(df_ct_tf, columns=["tf", "ct", "coef", "p",'SE', 'r_squared'])
#     df_ct_tf = df_ct_tf.sort_values(by=['tf', 'ct'])
#     df_ct_tf["p_adj"] = multipletests(df_ct_tf['p'], alpha=0.01, method="fdr_bh")[1]
#     df_ct_tf["negative_log_p_adj"] = -np.log10(df_ct_tf["p_adj"]+1e-10)
#     return df_ct_tf


def make_cor_dataframe(adata, adata_tfa, celltype_label='celltype'):
    A = adata_tfa.obsm[celltype_label].copy()
    A = (A-A.mean())/A.std()
    Y = adata_tfa.to_df()
    Y = Y-Y.mean()
    Y = Y/Y.std()
    mat_cor_tfa = Y.T.dot(A)/adata_tfa.n_obs
    
    tfs = np.intersect1d(adata_tfa.var_names, adata.var_names)
    Y = adata.to_df().loc[adata_tfa.obs_names, tfs]
    Y = Y-Y.mean()
    Y = Y/Y.std()
    mat_cor_rna = Y.T.dot(A)/adata_tfa.n_obs
    return mat_cor_tfa, mat_cor_rna
