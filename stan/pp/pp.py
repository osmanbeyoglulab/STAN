import pandas as pd
import numpy as np
# import NaiveDE
from scipy.linalg import lstsq, kron
from scipy.sparse.linalg import svds
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import jaccard
import networkx as nx
from matplotlib import pyplot as plt
import scanpy as sc
from sklearn.preprocessing import QuantileTransformer



def filter_tfs(tf_gene, min_jac=0.2, min_genes=7, min_tfs=10 ):
    #tf_gene=tf_gene.loc[tf_gene.index[tf_gene.T.sum()>min_genes],tf_gene.columns[tf_gene.sum()>min_tfs]]
    tf_gene.drop(tf_gene.index[tf_gene.T.sum()<min_genes], 0,inplace=True)
    tf_gene.drop(tf_gene.columns[tf_gene.sum()<min_tfs], 1,inplace=True)
    combine_tfs(tf_gene, cutoff=min_jac)

def combine_tfs(tf_gene, cutoff=0.5):
    '''
    aggrigates TFs that target mostly the same genes.
    Parameters
    ----------
    tf_gene
    cutoff

    Returns
    -------

    '''
    similar_tfs=tf_gene.corr(method=jaccard).melt(ignore_index=False).query("value<@cutoff")
    g = nx.Graph()
    for i in range(len(similar_tfs)):
        g.add_edge(similar_tfs.index.to_list()[i], similar_tfs.variable.to_list()[i])
    # add nodes/edges to graph

    d = list(nx.connected_components(g))
    # d contains disconnected subgraphs
    # d[0] contains the biggest subgraph

    for group in d:
        tf_gene['_'.join(list(group))]=(tf_gene.loc[:, list(group)].sum(axis=1)>1)+0
        tf_gene.drop(list(group), axis=1, inplace=True)

def filter_genes(adata, min_mean=0.5, min_var_mean_ratio=1, min_obs_ratio=0.2, filter_mt=True):
    '''
    Parameters
    ----------
    adata :
    min_mean :
    min_var_mean_ratio :
    min_obs_ratio :

    Returns
    -------
    '''

    sc.pp.filter_genes(adata, min_cells=round(min_obs_ratio*adata.n_obs), inplace=True)

    sc.pp.filter_genes(adata, min_counts=round(min_mean*adata.n_obs), inplace=True)
    counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)

    adata=adata[:,[a and b and c for a, b, c in zip((counts.mean()>min_mean).to_list(), (counts.mean()< min_var_mean_ratio*counts.std()**2).to_list(), ((counts>0).mean()>min_obs_ratio).to_list() )]]

    if filter_mt:
        adata=adata[:, ["MT-" not in x for x in adata.var_names]]
    return adata

def normalize_adata(adata, logscale=1e6, var_stable=True, quantile=True):

    counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
    adata.layers['raw_counts']=adata.X
    adata.layers["logp1"]=np.log(logscale*adata.layers["raw_counts"]/adata.layers["raw_counts"].sum(axis=1)+1)

    # if var_stable:
    #     sample_info = pd.DataFrame({'total_counts': counts.T.sum()})
    #     norm_expr = NaiveDE.stabilize(counts.T).T
    #     resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T.to_numpy()
    #     adata.layers['normalized']=resid_expr-resid_expr.mean(axis=1).reshape(-1,1)
    #     adata.layers['var_stablized']=norm_expr
    if quantile:
        qt = QuantileTransformer(output_distribution="normal")
        adata.layers['quantile_transformed']=qt.fit_transform(counts)


def common_genes(tf_gene, adata, min_cells=0):
    tf_to_keep=np.intersect1d(tf_gene.columns, adata.var_names)

    tf_to_keep=adata.var.loc[tf_to_keep].query("n_cells>@min_cells").index

    genes_to_keep=np.intersect1d(tf_gene.index, adata.var_names)

    tf_gene=tf_gene.loc[genes_to_keep, tf_to_keep]
    adata=adata[:,genes_to_keep]
    return (tf_gene, adata)

def assign_folds(adata, n_folds=5, train_percent=None):
    if train_percent is not None:
        adata.var["fold"]=np.random.uniform(0,1, adata.n_vars)<train_percent
    else:
        adata.var["fold"]=np.random.randint(0,n_folds, adata.n_vars)

import squidpy as sq
from sklearn.decomposition import PCA
from PIL import Image
# def make_kernel(adata,
#                 key=None,
#                 n=100,
#                 std=3,
#                 window_size=20,
#                 im_feats_weight=0.3,
#                 n_bins=20):
#
#     if key is None:
#         key=list(adata.uns['spatial'].keys())[0]
#
#     scale=adata.uns['spatial']["V1_Human_Lymph_Node"]['scalefactors']['tissue_hires_scalef']
#
#     image=Image.fromarray(np.uint8(adata.uns['spatial']["V1_Human_Lymph_Node"]["images"]['hires']*255))
#     image.convert("L")
#     image_feature_df=pd.DataFrame(None, columns=["mean"], index=adata.obs_names)
#
#     for i in (range(len(adata.obs_names))):
#         x=round(adata.obsm['spatial'][i,1]*scale)
#         y=round(adata.obsm['spatial'][i,0]*scale)
#         image_feature_df.iloc[i,0]=np.asarray(image.crop((y-window_size,x-window_size,y+window_size ,x+window_size))).mean()
#     adata.obsm['hist_feats']=image_feature_df.astype('float')
#     # sq_img = sq.im.ImageContainer(
#     #     adata.uns['spatial'][key]['images']['hires'],
#     #     scale=adata.uns['spatial'][key]['scalefactors']['tissue_hires_scalef']
#     # )
#     #
#     # sq.im.calculate_image_features(
#     #     adata,
#     #     sq_img,
#     #     features="histogram",
#     #     key_added="hist_feats",
#     #     show_progress_bar=False,
#     #     features_kwargs={"histogram": {"bins": n_bins}},
#     #     spot_scale=spot_scale
#     # )
#     #
#     # n_components = 2
#     # adata.obsm['hist_feats'] = pd.DataFrame(
#     #     PCA(n_components=n_components).fit_transform(adata.obsm['hist_feats']),
#     #     index=adata.obs_names,
#     #     columns=[str(i) for i in range(n_components)])
#
#     X = np.concatenate((adata.obsm['spatial'][:, 0:2], adata.obsm['hist_feats']), axis=1)
#     Xn = np.zeros(X.shape)
#     Xn[:, 0:2] = X[:, 0:2] - X[:, 0:2].mean()
#     Xn[:, 0:2] = Xn[:, 0:2] / np.max(np.abs(Xn[:, 0:2]))
#
#     if X.shape[1] > 2:
#         Xn[:, 2:X.shape[1]] = X[:, 2:X.shape[1]] - np.mean(X[:, 2:X.shape[1]], axis=0)
#         Xn[:, 2:X.shape[1]] = (im_feats_weight / np.sqrt(X.shape[1] - 2)) * np.divide(Xn[:, 2:X.shape[1]], np.std(
#             Xn[:, 2:X.shape[1]].astype(float), axis=0) + 1e-8)
#
#     omega = np.random.randn(X.shape[1], n) * std
#     proj = Xn.dot(omega)
#     adata.obsm['kernel'] = np.concatenate((np.cos(proj), np.sin(proj)), axis=1) / np.sqrt(n)

from scipy.spatial.distance import pdist, squareform
from PIL import Image
import squidpy as sq

def pixel_intensity(adata, key=None, window_size=20):
    if key is None:
        key=list(adata.uns['spatial'].keys())[0]

    scale=adata.uns['spatial'][key]['scalefactors']['tissue_hires_scalef']

    image=Image.fromarray(np.uint8(adata.uns['spatial'][key]["images"]['hires']*255))
    image.convert("L")
    adata.obsm['pixel_intensity']=pd.DataFrame(0, columns=["r", "g", "b"], index=adata.obs_names)
    adata.obs['pixel_intensity']=0
    for i in (range(len(adata.obs_names))):
        x=round(adata.obsm['spatial'][i,1]*scale)
        y=round(adata.obsm['spatial'][i,0]*scale)
        adata.obs['pixel_intensity'].iloc[i]=np.asarray(image.crop((y-window_size,x-window_size,y+window_size ,x+window_size))).mean()
        adata.obsm['pixel_intensity'].iloc[i, :]=np.asarray(image.crop((y-window_size,x-window_size,y+window_size ,x+window_size))).mean(axis=0).mean(axis=0).T
    adata.obsm['pixel_intensity']=(adata.obsm['pixel_intensity']-adata.obsm['pixel_intensity'].mean())/adata.obsm['pixel_intensity'].std()

def make_kernel(adata,
                n=100,
                banwidth=1,
                im_feats_weight=0.3):
    X = np.concatenate((adata.obsm['spatial'][:, 0:2], adata.obsm['pixel_intensity']), axis=1)

    x=X.T
    x=(x-x.mean(axis=1).reshape(-1,1))/x.std(axis=1).reshape(-1,1)
    Xn=x.T
    if X.shape[1] > 2:
        Xn[:, 2:X.shape[1]]=im_feats_weight*Xn[:, 2:X.shape[1]]
    pairwise_dists = squareform(pdist(Xn, 'euclidean'))
    adata.obsp['pairwise_dists']=pairwise_dists
    adata.obsp['kernel'] = (1/(np.sqrt(2)*np.pi*banwidth)**X.shape[1]) * np.exp(-pairwise_dists ** 2 /  (2 * banwidth ** 2))

    u,s,v=svds(adata.obsp['kernel'],n)
    adata.obsm['kernel']=u.dot(np.diag(s))

def make_kernel_RQ(adata,
                n=100,
                banwidth=1,
                im_feats_weight=0.3,
                   alpha=1):

    p=adata.obs['pixel_intensity'].to_numpy().reshape(-1,1)
    p=im_feats_weight*(p-p.mean())/p.std()

    X=adata.obsm['spatial']
    X=(X-X.mean())/X.std()

    X=np.concatenate((X, p), axis=1)
    pairwise_dists=pdist(X, 'sqeuclidean')
    pairwise_dists_mean=pairwise_dists.mean()
    pairwise_dists = squareform(pairwise_dists)/pairwise_dists_mean
    adata.obsp['pairwise_dists']=pairwise_dists

    tmp = (pairwise_dists) / (2 * alpha * banwidth**2)
    base = 1 + tmp
    adata.obsp['kernel'] = base**-alpha

    u,s,v=svds(adata.obsp['kernel'],n)
    adata.obsm['kernel']=u.dot(np.diag(s))


def kernel_smooth(X, K):
    K_row_norm=(K/K.sum(axis=0)).T
    return K_row_norm.dot(X)

def spatial_kernel_smooth(adata):
    adata.layers['spatial']=kernel_smooth(adata.to_df(), adata.obsp['kernel'])
    x=adata.layers['spatial'].T
    x=(x-x.mean(axis=1).reshape(-1,1))/x.std(axis=1).reshape(-1,1)
    adata.layers['spatial']=x.T