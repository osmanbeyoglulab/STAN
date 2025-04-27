import numpy as np
import pandas as pd
import time
from PIL import Image
from scipy.sparse.linalg import svds
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def recenter(X):
    x = X.T
    x = (x-x.mean(axis=1).reshape(-1,1))/x.std(axis=1).reshape(-1,1)
    Xn = x.T
    return Xn


def rbf_kernel(X):
    pw_dist = squareform(pdist(X, 'euclidean'))
    K = np.exp(-pw_dist**2 / X.shape[1])
    return K
    

def pixel_intensity(adata, key=None, windowsize=20):
    if key is None:
        key = list(adata.uns['spatial'].keys())[0]
    scale = adata.uns['spatial'][key]['scalefactors']['tissue_hires_scalef']

    image = np.uint8(adata.uns['spatial'][key]["images"]['hires']*255)
    if image.shape[-1]>3:
        image = image[:,:,:3]

    image = Image.fromarray(image)
    image.convert("L")
    intensity_3d = np.zeros((len(adata.obs_names), 3))
    # intensity_1d = np.zeros(len(adata.obs_names))
    x = np.round(adata.obsm['spatial'][:,1]*scale)
    y = np.round(adata.obsm['spatial'][:,0]*scale)
    for i in (range(len(adata.obs_names))):
        subimage = np.asarray(image.crop((y[i]-windowsize, x[i]-windowsize, y[i]+windowsize, x[i]+windowsize)))
        # intensity_1d[i] = subimage.mean()
        intensity_3d[i, :] = subimage.mean(axis=0).mean(axis=0).T
    # adata.obs['pixel'] = intensity_1d
    adata.obsm['pixel'] = ( intensity_3d - intensity_3d.mean(axis=0) ) / intensity_3d.std(axis=0)


def make_kernel(adata, X, n=100, kernel_name='kernel'):
    Xn = recenter(X)
    adata.obsp[kernel_name] = rbf_kernel(Xn)
    # pw_dist = squareform(pdist(Xn, 'euclidean'))
    # adata.obsp[kernel_name] = (1/(np.sqrt(2*np.pi)*bandwidth)**Xn.shape[1]) * np.exp(-pw_dist**2 / (2 * bandwidth**2))
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp[kernel_name], n)
    adata.obsm[kernel_name] = u.dot(np.diag(s))


def make_kernel_from_pixel(adata, n=100, kernel_name='kernel', im_feats_weight=0.3):
    Xn = recenter(np.concatenate((adata.obsm['spatial'][:, 0:2], adata.obsm['pixel']), axis=1))
    Xn[:, 2:Xn.shape[1]] = im_feats_weight*Xn[:, 2:Xn.shape[1]]
    adata.obsp[kernel_name] = rbf_kernel(Xn)
    # pw_dist = squareform(pdist(Xn, 'euclidean'))
    # adata.obsp[kernel_name] = (1/(np.sqrt(2*np.pi)*bandwidth)**Xn.shape[1]) * np.exp(-pw_dist**2 / (2 * bandwidth**2))
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp[kernel_name], n)
    adata.obsm[kernel_name] = u.dot(np.diag(s))


def get_feature(path, n_components=100, pct_variance=0.85):
    feature = np.load(path)
    feature = feature[:, ~np.all(feature == 0, axis=0)]
    n_components = 100
    pca = PCA(n_components=n_components)
    pca_feature = pca.fit_transform(feature)
    s = 0
    for i in range(n_components):
        s += pca.explained_variance_ratio_[i]
        if s>pct_variance:
            break
    n_features = i+1

    scaler = StandardScaler()
    norm_feature = scaler.fit_transform(pca_feature[:, :n_features])
    return norm_feature


def make_kernel_from_feature(adata, feature, kernel_name, n=250, im_feats_weight = 0.1):
    Xn = recenter(np.concatenate((adata.obsm['spatial'][:, 0:2], feature), axis=1))
    Xn[:, 2:] = im_feats_weight*Xn[:, 2:]
    adata.obsp[kernel_name] = rbf_kernel(Xn)
    # pw_dist = squareform(pdist(Xn, 'euclidean'))
    # adata.obsp[kernel_name] = (1/(np.sqrt(2*np.pi)*bandwidth)**Xn.shape[1]) * np.exp(-pw_dist**2 / (2 * bandwidth**2))
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp[kernel_name], n)
    adata.obsm[kernel_name] = u.dot(np.diag(s))
    return 