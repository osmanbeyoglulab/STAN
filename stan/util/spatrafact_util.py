import warnings

import pandas as pd
import numpy as np
from scipy.linalg import lstsq, kron
from scipy.sparse.linalg import svds
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import jaccard
import networkx as nx
from matplotlib import pyplot as plt
import scanpy as sc
from scipy.spatial import Voronoi
from scipy.stats import rankdata


def corr(A,B, axis=0, method="pearsonr"):
    if axis == 1:
        a=A.T
        b=B.T
    else:
        a=A
        b=B

    if method=="pearsonr":
        cor_method=pearsonr
    elif method=="spearmanr":
        cor_method=spearmanr
    else:
        print("Warning: unknown corr method, using pearsonr")

    if type(A) is pd.core.frame.DataFrame:
        return [pearsonr(B.iloc[:, i],A.iloc[:, i])[0] for i in range(A.shape[1])]
    else:
        return [pearsonr(B[:, i],A[:, i])[0] for i in range(A.shape[1])]






def adata2pd(adata, layer="normalized", im_feats=True):
    counts = adata.to_df(layer)
    if im_feats and 'hist_feats' in adata.obsm:
        X=np.concatenate((adata.obsm['spatial'][:,0:2], adata.obsm['hist_feats']), axis=1)
    else:
        X=adata.obsm['spatial'][:,0:2]
    X=X.astype('float32')
    return (counts, X)


def sci(A,B, X=None, W=None):

    if W is None:
        vor = Voronoi(X)
        W=vor.ridge_points
        W=np.vstack((W, np.stack((W[:,1], W[:,0])).T))

    elif W.shape[0]==W.shape[1]:
        W=np.stack(np.where(W), axis=0).T

    n_spots=np.max(W)+1

    if not (A.shape[0]==n_spots and B.shape[0]==n_spots):
        print("dimensions dont line up, dim 0 should be n_spots")


    if type(A) is pd.core.frame.DataFrame:
        x=A.to_numpy()
        y=B.to_numpy()
    else:
        x=A
        y=B

    # if method=="pearsonr":
    #     pass
    # elif method=="spearmanr":
    #     x=rankdata(x, axis=0)
    #     y=rankdata(y, axis=0)
    #
    # else:
    #     print("Warning: unknown corr method, using pearsonr")


    x_centered = x - np.mean(x, axis=0)
    y_centered = y - np.mean(y, axis=0)

    x_norm = np.sqrt(np.sum(x_centered ** 2, axis=0))
    y_norm = np.sqrt(np.sum(y_centered ** 2, axis=0))

    x_centered=x_centered/x_norm
    y_centered=y_centered/y_norm

    cor=x_centered[W[:,0],:].T.dot(y_centered[W[:,1],:])*n_spots/(2*len(W))

    if type(A) is pd.core.frame.DataFrame:
        cor=pd.DataFrame(cor, index=A.columns, columns=B.columns)

    return cor


def sci_test(A,B, X=None, W=None, n_trials=100):
    if type(A) is pd.core.frame.DataFrame:
        x=A.to_numpy()
        y=B.to_numpy()
    else:
        x=A
        y=B

    sci_true=sci(x,y, X=X, W=W)

    count=np.ones(sci_true.shape)
    for i in range(n_trials):
        count=count+(np.abs(sci_true)<=np.abs(sci(x[np.random.permutation(x.shape[0]),:],y[np.random.permutation(y.shape[0]),:], X=X, W=W)))

    return count/(n_trials+1)

def make_kernel(X, n, std=2, im_feats_weight=0.1):
    Xn=np.zeros(X.shape)
    Xn[:,0:2]=X[:,0:2]-X[:,0:2].mean()
    Xn[:,0:2]=Xn[:,0:2]/np.max(np.abs(Xn[:,0:2]))

    if X.shape[1]>2:
        Xn[:,2:X.shape[1]]=X[:,2:X.shape[1]]-np.mean(X[:,2:X.shape[1]], axis=0)
        Xn[:,2:X.shape[1]]=(im_feats_weight/np.sqrt(X.shape[1]-2))*np.divide(Xn[:,2:X.shape[1]],np.std(Xn[:,2:X.shape[1]].astype(float), axis=0)+1e-8)

    omega=np.random.randn(X.shape[1],n) * std
    proj=Xn.dot(omega)
    phi=np.concatenate((np.cos(proj), np.sin(proj)), axis=1)/np.sqrt(n)
    return phi


