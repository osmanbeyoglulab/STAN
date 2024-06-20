import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

figsize = 5
fontsize = 18


def plot_validation(adata):
    xstring = "pred_cor_ridge"
    ystring = "pred_cor_stan"

    plt.figure(figsize=(figsize, figsize), dpi=300)
    plt.rc('font', size=fontsize) 
    lim_min = np.minimum(adata.obs[xstring], adata.obs[ystring])
    lim_max = np.maximum(adata.obs[xstring], adata.obs[ystring])
    plt.plot([lim_min, lim_max], [lim_min, lim_max], '-', alpha=0.25, color='grey')
    sns.scatterplot(data=adata.obs, x=xstring, y=ystring, s=10, hue="n_counts", linewidth=0, palette='flare')
    plt.ylabel("STAN")
    plt.xlabel("Ridge")
    plt.legend(title="UMI Count", loc='upper right', bbox_to_anchor=(1.5, 1), columnspacing=0.5, ncol=1, handletextpad=0, frameon=False)
    plt.title('Cross Validation Performance\n(Pearson r)')


def plot_spatial_activity(adata, genes, points, edges):
    df = sc.get.rank_genes_groups_df(adata, group='GC')
    df.index = df['names']

    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(ngenes*figsize, figsize), dpi=300)
    plt.rc('font', size=fontsize) 
    for i in range(ngenes):
        tf = genes[i]
        sc.pl.spatial(adata, color=tf, 
            size=1.8, alpha_img=0, color_map="plasma", ax=axs[i], show=False, 
            legend_fontsize=fontsize, colorbar_loc='right')
        title = tf + ' activity\np_adj=%.2e'%df.loc[tf,'pvals_adj']
        axs[i].set_title(title, fontsize=fontsize)
        axs[i].set_xlabel("")
        axs[i].set_ylabel("")
        for ii, jj in edges:
            axs[i].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k', linewidth=0.9)
    plt.tight_layout(pad=0.6)


def plot_spatial_expression(adata, genes, points, edges):
    df = sc.get.rank_genes_groups_df(adata, group='GC')
    df.index = df['names']

    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(ngenes*figsize, figsize), dpi=300)
    plt.rc('font', size=fontsize) 
    for i in range(ngenes):
        tf = genes[i]
        sc.pl.spatial(adata, color=tf, 
            size=1.8, alpha_img=0, color_map="viridis", ax=axs[i], show=False, 
            legend_fontsize=fontsize, colorbar_loc='right')
        title = tf + ' mRNA expr\np_adj=%.2e'%df.loc[tf,'pvals_adj']
        axs[i].set_title(title, fontsize=fontsize)
        axs[i].set_xlabel("")
        axs[i].set_ylabel("")
        for ii, jj in edges:
            axs[i].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'w', linewidth=0.9)
    plt.tight_layout(pad=0.6)


def plot_umap(adata, palette=None, is_tf=True):
    if is_tf:
        title = 'Clustering based on TF Activity'
    else: 
        title = 'Clustering based on mRNA Expression'
    fig, axs = plt.subplots(1,2, figsize=(figsize*2.25, figsize), width_ratios=[1.25,1], dpi=300)
    plt.rc('font', size=fontsize) 
    sc.pl.umap(adata, color="leiden", size=40, palette=palette, ax=axs[0], 
               show=False, frameon=True)
    sc.pl.spatial(adata, color="leiden", size=1.8, alpha_img=0, palette=palette, ax=axs[1], 
                  show=False, frameon=True, legend_fontsize=fontsize)
    for ax in axs.flatten():
        ax.set_ylabel(ax.get_ylabel(),labelpad=-1)
        ax.set_xlabel(ax.get_xlabel(),labelpad=-1)
        ax.set_title("")
    axs[0].legend().remove()
    axs[1].legend(title='Leiden\nCluster', loc='upper right', bbox_to_anchor=(1.5, 1), columnspacing=0.5, ncol=2, handletextpad=0, frameon=False)
    plt.suptitle(title, fontsize=fontsize*1.2)


# ========== ========== ========== ========== ========== ========== ========== ========== ========== ==========


def make_ct_tf_dataframe(adata, celltype_label='celltype'):
    df_ct_tf = []
    A = adata.obsm[celltype_label].copy()
    Y = adata.to_df()
    Y = Y-Y.mean()
    Y = Y/Y.std()
    for tf in adata.var_names:
        y = Y[tf]
        results = sm.OLS(y,A).fit()
        for ct in A.columns:
            p = results.pvalues[ct]
            coef = results.params[ct]
            se = results.bse[ct]
            df_ct_tf.append([tf, ct, coef, p, se, results.rsquared])

    df_ct_tf = pd.DataFrame(df_ct_tf, columns=["tf", "ct", "coef", "p",'SE', 'r_squared'])
    df_ct_tf = df_ct_tf.sort_values(by=['tf', 'ct'])
    df_ct_tf["p_adj"] = multipletests(df_ct_tf['p'], alpha=0.01, method="fdr_bh")[1]
    df_ct_tf["negative_log_p_adj"] = -np.log10(df_ct_tf["p_adj"]+1e-10)
    return df_ct_tf


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


def plot_heatmap(df_ct_tf, tf_list, ct_list, clip=10):
    data = df_ct_tf.query("tf in @tf_list and ct in @ct_list")
    data.columns = ['tf', 'ct', 'TF Score', 'p', 'SE', 'r_squared', 'p_adj', '-log(p_adj)']
    x = 'tf'
    y = 'ct'

    data[x] = data[x].astype("category")
    data[y] = data[y].astype("category")
    x_lab = data[x].cat.categories
    y_lab = data[y].cat.categories

    f = sns.clustermap(data.pivot(index=y, columns=x, values="TF Score"),figsize=(0.1,0.1), cmap='PiYG')
    x_lab = x_lab[f.dendrogram_col.reordered_ind]
    y_lab = y_lab[f.dendrogram_row.reordered_ind]

    data[x] = data[x].cat.reorder_categories(x_lab)
    data[y] = data[y].cat.reorder_categories(y_lab)
    data = data.sort_values([x, y])
    data["TF Score"] = data["TF Score"].clip(-clip, clip)

    figsize = 0.25
    plt.figure(figsize=(figsize*len(x_lab)*0.9, figsize*len(y_lab)), dpi=300)
    plt.rc('font', size=fontsize/1.8) 
    ax = sns.scatterplot(data=data,x=x, y=y, palette="PiYG_r", hue="TF Score", size="-log(p_adj)")
    plt.legend(bbox_to_anchor=(1.2,1.), loc='upper right', 
        columnspacing=0.5, handletextpad=0, frameon=False, fontsize=fontsize/1.8)

    ax.set_xticklabels(x_lab, rotation = 90)
    ax.set_xlim(-0.5, -0.5+len(x_lab))
    ax.set_ylim(-0.5, -0.5+len(y_lab))
    ax.set_xlabel("")
    ax.set_ylabel("")


# ========== ========== ========== ========== ========== ========== ========== ========== ========== ==========


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


def plot_heatmap_correlation(data, title):
    n = len(data['sample'].unique())-1
    plt.figure(figsize=(figsize,figsize), dpi=300)
    plt.rc('font', size=fontsize)
    ax = sns.heatmap(data.pivot(columns="sample", index="TF", values="TFa").dropna().corr().round(2).iloc[1:,:-1],
                cmap="coolwarm", center=0, annot=True, mask=np.tril(np.ones(n)-np.eye(n)).T, vmax=1, vmin=-1)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0, ha='right')
    plt.title(title)
    plt.legend(title='Pearson r', frameon=False)
