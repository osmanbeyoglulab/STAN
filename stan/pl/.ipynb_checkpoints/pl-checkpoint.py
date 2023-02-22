import numpy as np
from matplotlib import pyplot as plt
from spatrafact.util.spatrafact_util import make_kernel, adata2pd
from scipy.stats import spearmanr
import pandas as pd
from sklearn.decomposition import FastICA
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering

def plot_kernel_samples(X, n_samples, std=2, im_feat_weights=0.1):
    kernel=make_kernel(X, n_samples, std=std, im_feats_weights=im_feat_weights)
    x=np.random.randn(8,kernel.shape[1]).dot(kernel.T)
    for i in range(8):
        plt.subplot(2,4,i+1)
        plt.scatter(X[:,0], -X[:,1], c=x[i,:], s=10)


def plot_hist_feats(adata, n_feats=None):
    _,X=adata2pd(adata)
    if n_feats is None:
        n_feats=adata.obsm['hist_feats'].shape[1]
    for i in range(n_feats):
        plt.subplot(np.ceil(n_feats/2), 4, i+1)
        plt.scatter(X[:,0],-X[:, 1], c=X[:, 2+i], s=6)

def plot_gene_cor(adata, layer1, layer2, gene):
    x1=adata.obs_vector(gene, layer=layer1)
    x2=adata.obs_vector(gene, layer=layer2)
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    plt.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1], c=x1, s=3)
    plt.title(gene+" "+layer1)
    plt.subplot(1,3,2)
    plt.title(gene+" "+layer2)
    plt.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1], c=x2, s=3)
    plt.subplot(1,3,3)
    plt.scatter(x1,x2, s=0.1)
    plt.xlabel(gene+" "+layer1)
    plt.ylabel(gene+" "+layer2)
    plt.title("Spearman Corr: "+str(round(spearmanr(x1,x2)[0], 3)))

def spatial_patterns(adata, layer='tf_activities', n_components=9, plot_results=True, s=3, n_tfs_to_plot=5, n_tfs=10):
    W=adata.obsm[layer].T

    fica=FastICA(n_components=n_components)
    patterns=fica.fit_transform(W.T)

    adata.obsm["patterns"]=patterns
    components=np.divide(fica.components_, np.linalg.norm(fica.components_, axis=0))

    components=pd.DataFrame(data=components[:,0:-1], columns=adata.uns["tf_names"])
    adata.uns["components"]=components

    top_tfs=pd.DataFrame([components.T.sort_values(by=i, ascending=False).index[0:n_tfs].to_list() for i in range(n_components)], index=["Factor "+str(i) for i in range(n_components)])

    if plot_results:
        plt.figure(figsize=(4*np.ceil(n_components/4), 8))
        for i in range(0,n_components):
            plt.subplot(np.ceil(n_components/4),4,i+1)
            plt.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1],  c=patterns[:,i], s=1)
            plt.title("Factor "+str(i))
            plt.gca().get_xaxis().set_ticks([])
            plt.gca().get_yaxis().set_ticks([])
        plt.suptitle("Spatial Patterns")

        plt.figure(figsize=( (1+n_tfs_to_plot*2), n_components*2))
        for i in range(0,n_components):
            plt.subplot(n_components, 1+n_tfs_to_plot, (n_tfs_to_plot+1)*i+1)
            plt.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1],  c=patterns[:,i], s=0.5)
            plt.title("Factor "+str(i))
            plt.gca().get_xaxis().set_ticks([])
            plt.gca().get_yaxis().set_ticks([])

            for j in range(n_tfs_to_plot):
                plt.subplot(n_components, 1+n_tfs_to_plot, (n_tfs_to_plot+1)*i+1+j+1)
                plt.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1],  c=W.loc[:,top_tfs.iloc[i,j]], s=0.5)
                plt.title(top_tfs.iloc[i,j])
                plt.gca().get_xaxis().set_ticks([])
                plt.gca().get_yaxis().set_ticks([])
        plt.suptitle("Top TFs for each pattern")
        print("Top TFs in each pattern:")
        print(top_tfs.to_string(header=False))

def spatial_cluster(adata,layer='tf_actvities', n_clusters=5, plot_results=True, s=3,figsize=(10,5)):
    W=adata.obsm[layer]
    #I tried adding the coordinates as a feature for clustering but Hatice said we probably shouldnt add that as a bias because this is just for evaluating
    #kmeans = SpectralClustering(n_clusters=6, random_state=0).fit(np.concatenate((W.T,spatial_weight* self.X[:,0:2]/np.max(self.X[:,0:2])), axis=1))
    kmeans = SpectralClustering(n_clusters=n_clusters, random_state=0).fit(W)
    clusters=kmeans.labels_
    W_embedded = TSNE(n_components=2,init='random').fit_transform(W.T)

    if plot_results:
        fig, (ax1, ax2, ax3)=plt.subplots(1,3, figsize=figsize)
        ax1.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1], c=clusters, s=s)
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        ax1.set_ylabel("x coordinate")
        ax1.set_xlabel("y coordinate")

        tsne_scatter=ax2.scatter(W_embedded[:,0], W_embedded[:,1], c=clusters, s=s)
        ax2.get_xaxis().set_ticks([])
        ax2.get_yaxis().set_ticks([])
        ax2.set_ylabel("TSNE2")
        ax2.set_xlabel("TSNE1")

        ax3.legend(*tsne_scatter.legend_elements(), bbox_to_anchor=(0.5, 1), title="Cluster")
        ax3.axis('off')

        ax1.set_box_aspect(1)
        ax2.set_box_aspect(1)
        ax3.set_box_aspect(1)


def spatial(adata, vars, layer= ['raw_counts'], n_cols=4, s=3):
    plt.figure(figsize=(20, 5*np.ceil(len(vars) / n_cols)))
    for i, var in enumerate(vars):
        plt.subplot(np.ceil(len(vars) / n_cols), 4, i+1)
        plt.scatter(adata.obsm['spatial'][:,0], -adata.obsm['spatial'][:,1], c=adata.obs_vector(var,layer), s=s)
        plt.title(var)
