import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import Delaunay


figsize = 2.5
fontsize = 9
dpi = 150


def plot_heatmap_correlation(data, title):
    n = len(data['sample'].unique())-1
    plt.figure(figsize=(figsize, figsize), dpi=dpi)
    plt.rc('font', size=fontsize)
    ax = sns.heatmap(data.pivot(columns="sample", index="TF", values="TFa").dropna().corr().round(2).iloc[1:,:-1],
                cmap="coolwarm", center=0, annot=True, mask=np.tril(np.ones(n)-np.eye(n)).T, vmax=1, vmin=-1, fmt='.2f')
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0, ha='right')
    plt.title(title)
    plt.legend(title='Pearson r', frameon=False)


def make_pat_tf_dataframe(adata_list):
    dfs_tf_pat = dict()
    for sample in adata_list.keys():
        adata = adata_list[sample]
        categories = adata.obs['pathology'].cat.categories
        mean_df = adata.to_df().groupby(adata.obs['pathology']).mean()

        df = []
        for i in range(len(adata.uns['rank_genes_groups']['names'])):
            for j, cat in enumerate(categories):
                tf = adata.uns['rank_genes_groups']['names'][i][j]
                pval = adata.uns['rank_genes_groups']['pvals_adj'][i][j]
                mean = mean_df.loc[cat, tf]
                df.append([sample, tf, cat, mean, -np.log10(pval+1e-10), i])
                
        dfs_tf_pat[sample] = pd.DataFrame(df, columns=['sample', 'TF', 'Pathology', 'TFa', '-log(p_adj)', 'rank'])
        dfs_tf_pat[sample]['abs_TFa'] = dfs_tf_pat[sample]['TFa'].abs()
    return dfs_tf_pat


def find_edges(adata, sample, category, r=10, alpha=20, only_outer=True):
    points = np.asarray(adata[adata.obs['pathology']==category].obsm['spatial']*adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'])
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


def dendrogram_order(data):
    x = "TF"
    y = 'Pathology'
    
    data[x]=data[x].astype("category")
    data[y]=data[y].astype("category")
    x_lab=data[x].cat.categories
    y_lab=data[y].cat.categories

    f = sns.clustermap(data.pivot(index=y, columns=x, values="TFa"),figsize=(0.1,0.1), cmap='PiYG')
    x_lab = x_lab[f.dendrogram_col.reordered_ind]
    y_lab = y_lab[f.dendrogram_row.reordered_ind]
    print(x_lab)
    print(y_lab)
    plt.close()
    return