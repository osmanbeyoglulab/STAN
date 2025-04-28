import numpy as np
import pandas as pd
import time
from PIL import Image
from scipy.sparse.linalg import svds
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def recenter(X):
    """
    Standardize the input data matrix by centering and scaling features (columns).
    
    For each feature (column) in the input matrix, this function:
    1. Centers the data by subtracting the mean (mean normalization)
    2. Scales the data by dividing by the standard deviation (unit variance)
    
    Parameters:
    -----------
    X : numpy.ndarray or array-like
        Input data matrix with shape (n_samples, n_features)
        
    Returns:
    --------
    numpy.ndarray
        Standardized data matrix with zero mean and unit variance for each feature
    """
    x = X.T
    x = (x-x.mean(axis=1).reshape(-1,1))/x.std(axis=1).reshape(-1,1)
    Xn = x.T
    return Xn


def rbf_kernel(X):
    """
    Compute the Radial Basis Function (RBF) kernel matrix for the input data.
    """
    pw_dist = squareform(pdist(X, 'euclidean'))
    K = np.exp(-pw_dist**2 / X.shape[1])
    return K
    

def pixel_intensity(adata, key=None, windowsize=20):
    """
    Compute normalized pixel intensity features around each spot in spatial transcriptomics data.
    
    For each spot, extracts a window around its location in the H&E image, computes average RGB intensities,
    and returns standardized intensity features (z-scored per channel).

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing spatial coordinates and image data
    key : str, optional
        Key identifying which tissue sample to use (default: first available sample)
    windowsize : int
        Half-width of the square window to extract around each spot (in pixels)
        Total window size will be (2*windowsize+1) x (2*windowsize+1)

    Returns:
    --------
    None (modifies adata in-place by adding 'pixel' to obsm)
    """
    
    if key is None:
        key = list(adata.uns['spatial'].keys())[0]

    # Get scaling factor to convert spot coordinates to image pixels
    scale = adata.uns['spatial'][key]['scalefactors']['tissue_hires_scalef']

    # Extract high-resolution image and convert to 8-bit RGB
    image = np.uint8(adata.uns['spatial'][key]["images"]['hires']*255)

    # Ensure image has only 3 channels (drop alpha channel if present)
    if image.shape[-1]>3:
        image = image[:,:,:3]

    # Convert to PIL Image object for cropping operations
    image = Image.fromarray(image)
    image.convert("L")

    # Initialize array to store RGB intensities for each spot
    intensity_3d = np.zeros((len(adata.obs_names), 3))

    # Get spot coordinates scaled to image resolution
    x = np.round(adata.obsm['spatial'][:,1]*scale)
    y = np.round(adata.obsm['spatial'][:,0]*scale)
    for i in (range(len(adata.obs_names))):
        subimage = np.asarray(image.crop((y[i]-windowsize, x[i]-windowsize, y[i]+windowsize, x[i]+windowsize)))
        intensity_3d[i, :] = subimage.mean(axis=0).mean(axis=0).T
    adata.obsm['pixel'] = ( intensity_3d - intensity_3d.mean(axis=0) ) / intensity_3d.std(axis=0)


def make_kernel(adata, X, n=100, kernel_name='kernel'):
    """
    Create and store a kernel matrix and its SVD decomposition in an AnnData object.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object to store the results in.
    X : array-like
        Input data matrix (cells x features).
    n : int, optional (default=100)
        Number of singular vectors/values to compute. Will be reduced if larger than number of samples.
    kernel_name : str, optional (default='kernel')
        Base name to use for storing results in adata (will use obsp and obsm).
    """
    Xn = recenter(X)
    adata.obsp[kernel_name] = rbf_kernel(Xn)
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp[kernel_name], n)
    adata.obsm[kernel_name] = u.dot(np.diag(s))


def make_kernel_from_pixel(adata, n=100, kernel_name='kernel', im_feats_weight=0.3):
    Xn = recenter(np.concatenate((adata.obsm['spatial'][:, 0:2], adata.obsm['pixel']), axis=1))
    Xn[:, 2:Xn.shape[1]] = im_feats_weight*Xn[:, 2:Xn.shape[1]]
    adata.obsp[kernel_name] = rbf_kernel(Xn)
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp[kernel_name], n)
    adata.obsm[kernel_name] = u.dot(np.diag(s))


def get_feature(path, n_components=100, pct_variance=0.85):
    """
    Load feature data, preprocess it, perform PCA, and return normalized components.
    
    Parameters:
    -----------
    path : str
        Path to the numpy file containing the feature data.
    n_components : int, optional (default=100)
        Maximum number of PCA components to compute.
    pct_variance : float, optional (default=0.85)
        Target percentage of variance to retain with PCA components.
        
    Returns:
    --------
    numpy.ndarray
        Normalized feature matrix with selected PCA components.
    """
    
    # Load feature data from numpy file
    feature = np.load(path)
    
    # Remove columns (features) that are all zeros
    feature = feature[:, ~np.all(feature == 0, axis=0)]
    
    # Initialize PCA with specified number of components
    pca = PCA(n_components=n_components)
    
    # Fit PCA and transform the data
    pca_feature = pca.fit_transform(feature)
    
    # Calculate cumulative explained variance to find the number of components
    # needed to reach the target variance percentage
    s = 0  # Cumulative variance tracker
    for i in range(n_components):
        s += pca.explained_variance_ratio_[i]
        if s > pct_variance:
            break  # Stop when we reach target variance
    
    # The actual number of features we'll keep (+1 because of 0-based indexing)
    n_features = i + 1
    
    # Standardize the selected PCA components (mean=0, variance=1)
    scaler = StandardScaler()
    norm_feature = scaler.fit_transform(pca_feature[:, :n_features])
    
    return norm_feature


def make_kernel_from_feature(adata, feature, kernel_name, n=250, im_feats_weight = 0.1):
    Xn = recenter(np.concatenate((adata.obsm['spatial'][:, 0:2], feature), axis=1))
    Xn[:, 2:] = im_feats_weight*Xn[:, 2:]
    adata.obsp[kernel_name] = rbf_kernel(Xn)
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp[kernel_name], n)
    adata.obsm[kernel_name] = u.dot(np.diag(s))
    return 
