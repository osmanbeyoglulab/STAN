import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import sklearn
import squidpy as sq
import statsmodels.api as sm
from anndata import AnnData
from path import Path
from sklearn.preprocessing import OrdinalEncoder
from statsmodels.stats.multitest import multipletests

figsize = 2.5
fontsize = 9
dpi = 150


def t_stat(adata, r_key):
    """
    Calculates t-statistics from correlation coefficients stored in AnnData object.
    
    Converts correlation values (r) to t-statistics using the formula:
    t = r * sqrt((n-2)/(1-r²))
    
    Args:
        adata (AnnData): Annotated data object containing:
            - varm['gene_tf']: Gene-TF interaction matrix (features x TFs)
            - obsm[r_key]: Correlation coefficients matrix (spots x TFs)
        r_key (str): Key in adata.obsm containing correlation coefficients
    
    Returns:
        numpy.ndarray: Matrix of t-statistics with same shape as input correlations
    """
    
    # Get gene-TF interaction matrix (features x transcription factors)
    A = adata.varm['gene_tf']
    
    # Get expression matrix transposed (genes x spots)
    b = adata.to_df().T
    
    # Extract correlation coefficients (spots x TFs)
    r = adata.obsm[r_key]
    
    # Calculate degrees of freedom (n_features - 2)
    # Where n_features = number of genes (from gene-TF matrix)
    n_samples = b.shape[1]    # Number of spots
    n_features, n_fsets = A.shape  # Genes x TFs dimensions
    df = n_features - 2       # Degrees of freedom for t-test
    
    # Convert correlations to t-statistics:
    # t = r * sqrt(df / (1 - r²))
    # Added small epsilon (1e-16) for numerical stability near r=±1
    x = r * np.sqrt(df / ((1.0 - r + 1.0e-16) * (1.0 + r + 1.0e-16)))
    
    return x


def annotate_lymphnode(adata, fpath="resources/lymphnode_annotation"):
    """
    Annotates lymph node spatial transcriptomics data with cell type proportions 
    and germinal center (GC) identification.
    
    Args:
        adata (AnnData): Spatial transcriptomics data object
        fpath (str): Path to directory containing annotation files
        
    Returns:
        AnnData: Modified AnnData object with added annotations:
            - obsm['celltype']: Normalized cell type proportions
            - obsm['celltype_raw']: Raw cell type densities
            - obs['germinal_center']: GC classification ("GC" or "Other")
    """
    
    # 1. Load Annotation Files
    # -----------------------
    # Read cell type density estimates (spot x cell type matrix)
    celltypes = pd.read_csv(Path(fpath) / "W_cell_density.csv", index_col=0)
    # Read manual germinal center annotations (binary: 1=GC, 0=non-GC)
    gc_annotation = pd.read_csv(Path(fpath) / "manual_GC_annot.csv", 
                              index_col=0).fillna(0).replace("GC", 1)

    # 2. Align Data with Annotations
    # -----------------------------
    # Find overlapping spots between data and annotations
    obs_names = np.intersect1d(celltypes.index, adata.obs_names)
    # Subset both objects to matching spots only
    adata = adata[obs_names]
    celltypes = celltypes.loc[obs_names]
    gc_annotation = gc_annotation.loc[obs_names]

    # 3. Store Cell Type Information
    # ----------------------------
    # Add raw cell type densities to obsm
    adata.obsm['celltype'] = celltypes
    # Clean up column names by removing prefix
    adata.obsm['celltype'].columns = [
        x.replace('mean_spot_factors','') 
        for x in adata.obsm['celltype'].columns
    ]
    # Create backup of raw values before normalization
    adata.obsm['celltype_raw'] = adata.obsm['celltype'].copy()
    # Normalize to proportions (sum to 1 per spot)
    adata.obsm['celltype'] = adata.obsm['celltype'].divide(
        adata.obsm['celltype'].sum(axis=1), 
        axis=0
    )

    # 4. Add Germinal Center Annotations
    # --------------------------------
    # Store binary GC annotations
    adata.obs['germinal_center'] = gc_annotation
    # Convert numeric labels to descriptive strings
    adata.obs['germinal_center'] = adata.obs['germinal_center'].map(
        {0: "Other", 1: "GC"}
    )

    return adata


def plot_validation(adata, xstring="pred_cor_ridge", ystring="pred_cor_stan", 
                    xlabel='Ridge', ylabel='STAN', title='Cross Validation Performance\n(Pearson r)'):
    plt.figure(figsize=(figsize, figsize), dpi=dpi)
    plt.rc('font', size=fontsize) 
    lim_min = np.minimum(adata.obs[xstring], adata.obs[ystring])
    lim_max = np.maximum(adata.obs[xstring], adata.obs[ystring])
    plt.plot([lim_min, lim_max], [lim_min, lim_max], '-', alpha=0.25, color='grey')
    sns.scatterplot(data=adata.obs, x=xstring, y=ystring, s=5, hue="n_counts", linewidth=0, palette='flare')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(title="UMI Count", loc='upper right', bbox_to_anchor=(1.4, 1), columnspacing=0.5, ncol=1, handletextpad=0, frameon=False)
    plt.title(title)


def plot_umap(adata, palette=None, title='Clustering'):
    fig, axes = plt.subplots(1,2, figsize=(figsize*2, figsize), dpi=dpi)

    plt.rc('font', size=fontsize) 
    sc.pl.umap(adata, color="leiden", size=20, palette=palette, ax=axes[0], 
               show=False, frameon=True)
    sc.pl.spatial(adata, color="leiden", size=1.8, alpha_img=0, palette=palette, ax=axes[1], 
                  show=False, frameon=True, legend_fontsize=fontsize)
    for ax in axes.flatten():
        ax.set_ylabel(ax.get_ylabel(),labelpad=-1)
        ax.set_xlabel(ax.get_xlabel(),labelpad=-1)
        ax.set_title("")
    axes[0].legend().remove()
    axes[1].legend(title='Leiden\nCluster', loc='upper right', bbox_to_anchor=(1.4, 1), columnspacing=0.5, ncol=1, handletextpad=0, frameon=False)
    plt.tight_layout(pad=1.5)
    plt.suptitle(title, fontsize=fontsize*1.2)


def get_activity(adata, key='tfa_stan'):
    s = key[:2]
    adata_tfa = AnnData(
        X = adata.obsm[key],
        obs = adata.obs,
        obsm = {name: obj for (name, obj) in adata.obsm.items() if s not in name},
        layers = {name: obj for (name, obj) in adata.obsm.items() if s in name})
    adata_tfa.uns = adata.uns
    return adata_tfa


def compute_spatial_expression(adata, genes):
    """
    Compute spatially weighted expression profiles for given genes.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix containing spatial gene expression data
    genes : list
        List of gene names to compute spatial expression for
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing spatially smoothed expression values for the specified genes
    """
    sq.gr.spatial_neighbors(adata, n_rings=1)
    A = adata.obsp['spatial_connectivities']
    mat = sc.get.obs_df(adata, genes)
    mat_neighbor = pd.DataFrame(
        (A + np.eye(adata.n_obs)).dot(mat)/(1+A.sum(axis=1)),
        index = mat.index,
        columns = mat.columns)
    return mat_neighbor


def compute_ari(adata, cluster_1, cluster_2):
    """
    Compute the Adjusted Rand Index (ARI) between two clustering results.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix containing cluster labels
    cluster_1/cluster_2 : str
        Column name in adata.obs for the first/secent clustering result
        
    Returns:
    --------
    float
        Adjusted Rand Index score between the two clusterings
    """
    label_1 = np.array(adata.obs[cluster_1]).reshape(-1,1)
    sklearn_encoder = OrdinalEncoder()
    encoder_1 = sklearn_encoder.fit_transform(label_1)

    label_2 = np.array(adata.obs[cluster_2]).reshape(-1,1)
    sklearn_encoder = OrdinalEncoder()
    encoder_2 = sklearn_encoder.fit_transform(label_2)

    s = sklearn.metrics.adjusted_rand_score(encoder_1.flatten(), encoder_2.flatten())
    return s


def plot_spatial_activity(adata, genes, points, edges):
    df = sc.get.rank_genes_groups_df(adata, group='GC')
    df.index = df['names']

    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(ngenes*figsize, figsize), dpi=dpi)
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
            axs[i].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k', linewidth=0.8)
    plt.tight_layout(pad=0.6)


def plot_spatial_expression(adata, genes, points, edges):
    df = sc.get.rank_genes_groups_df(adata, group='GC')
    df.index = df['names']

    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(ngenes*figsize, figsize), dpi=dpi)
    plt.rc('font', size=fontsize) 
    for i in range(ngenes):
        tf = genes[i]
        sc.pl.spatial(adata, color=tf, 
            size=1.8, alpha_img=0, color_map="viridis", ax=axs[i], show=False, 
            legend_fontsize=fontsize, colorbar_loc='right')
        title = tf + ' mRNA expression\np_adj=%.2e'%df.loc[tf,'pvals_adj']
        axs[i].set_title(title, fontsize=fontsize)
        axs[i].set_xlabel("")
        axs[i].set_ylabel("")
        for ii, jj in edges:
            axs[i].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'w', linewidth=0.8)
    plt.tight_layout(pad=0.6)
