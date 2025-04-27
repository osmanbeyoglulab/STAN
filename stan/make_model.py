import numpy as np
import pandas as pd
import scanpy as sc
import time
from scipy.stats import spearmanr, pearsonr

def assign_folds(adata, n_folds=5, train_percent=None, random_seed=0):
    """
    Assign genes to folds for cross-validation.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    n_folds : int
        Number of cross-validation folds (ignored if train_percent specified)
    train_percent : float, optional
        If specified, creates a single train/test split with this percentage as train
    random_seed : int
        Seed for reproducible random assignments
    """
    np.random.seed(random_seed)
    if train_percent is not None:
        # Binary split: True for train genes, False for test genes
        adata.var["fold"] = np.random.uniform(0, 1, adata.n_vars) < train_percent
    else:
        # Assign each gene to one of n_folds
        adata.var["fold"] = np.random.randint(0, n_folds, adata.n_vars)

class ModelBase:
    """Base class for spatial gene expression modeling"""
    
    def __init__(self, adata, kernel_name='kernel', layer="dca", gene_tf=None, intercept=False):
        """
        Initialize base model with data and parameters.
        
        Parameters:
        -----------
        adata : AnnData
            Spatial transcriptomics data
        kernel_name : str
            Key in adata.obsm containing spatial kernel matrix
        layer : str
            Layer in adata containing expression data
        gene_tf : str, optional
            Key in adata.varm containing gene-TF relationships
        intercept : bool
            Whether to include intercept term
        """
        self.adata = adata
        self.intercept = intercept

        # Load gene-TF relationships (try alternative name if primary not found)
        try:
            gene_tf = adata.varm['gene_tf']
        except:
            gene_tf = adata.varm['gene_pw']

        self.D = gene_tf  # Gene-TF matrix
        self.Y = adata.to_df(layer).T  # Expression matrix (genes x spots)
        
        # Align genes between expression and TF data
        self.genes = np.intersect1d(self.Y.index, self.D.index)
        self.D = self.D.loc[self.genes]
        self.Y = self.Y.loc[self.genes]

        # Set up cross-validation folds
        self.n_folds = np.max(adata.var['fold']) + 1
        self.n_genes = len(self.genes)
        self.n_spots = adata.n_obs
        self.n_tfs = self.D.shape[1]

        # Create train/test splits for each fold
        self.training_genes = [adata.var.query("fold != @i").index for i in range(self.n_folds)]
        self.testing_genes = [adata.var.query("fold == @i").index for i in range(self.n_folds)]
        self.D_train = [self.D.loc[self.training_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.D_test = [self.D.loc[self.testing_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.Y_train = [self.Y.loc[self.training_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.Y_test = [self.Y.loc[self.testing_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.D = self.D.to_numpy()  # Full gene-TF matrix as numpy array

        # Load spatial kernel matrix
        self.K = adata.obsm.get(kernel_name)
        
        # Precompute SVD decompositions for each training set
        self.svdD = []
        for i in range(self.n_folds):
            self.svdD.append(np.linalg.svd(self.D_train[i], full_matrices=False))
        self.svdK = None  # Will be computed for spatial models

        # Initialize model weights and predictions
        self.W = [None] * self.n_folds
        self.y_pred = [None] * self.n_folds

    def evaluate(self, fold=0, gene_set="testing"):
        """
        Evaluate model performance for a given fold.
        
        Parameters:
        -----------
        fold : int
            Fold to evaluate (-1 for concatenated predictions)
        gene_set : str
            "testing" or "training" to evaluate on test/train genes
            
        Returns:
        --------
        tuple: (spot_correlations, gene_correlations)
        """
        if fold == -1:
            # Use concatenated predictions if fold is -1
            y_pred = np.asarray(self.y_pred_concat)
            Y = np.asarray(self.Y)
        else:
            # Get weights and appropriate data split
            W = self.W[fold]
            if gene_set == "testing":
                Y = self.Y_test[fold]
                D = self.D_test[fold]
            elif gene_set == "training":
                Y = self.Y_train[fold]
                D = self.D_train[fold]
            y_pred = D.dot(W)  # Compute predictions

        # Compute Pearson correlations
        y_pred = y_pred.astype(float)
        Y = Y.astype(float)

        # Correlation per spot and per gene
        cor = [pearsonr(y_pred[:, spot], Y[:, spot])[0] for spot in range(self.n_spots)]
        gene_cor = [pearsonr(y_pred[i, :].T, Y[i, :].T)[0] for i in range(Y.shape[0])]
        return cor, gene_cor

    def fit(self, grid_search_params=None, fixed_params=dict(), verbose=False, stages=5, n_steps=5, axis=0):
        """Fit model with optional hyperparameter tuning."""
        t1 = time.time()
        if grid_search_params is None:
            # Use fixed parameters if no grid search specified
            self.params = fixed_params
        else:
            # Perform grid search for parameter tuning
            self.params = self.grid_search(n_steps, grid_search_params, 
                                         fixed_params=fixed_params, 
                                         verbose=verbose, 
                                         stages=stages, 
                                         axis=axis)
        
        # Train model on each fold with best parameters
        params = self.params
        for fold in range(self.n_folds):
            self.train(fold=fold, **params)
            self.y_pred[fold] = self.D_test[fold].dot(self.W[fold])

        # Combine predictions across folds
        self.W_concat = np.mean(self.W, axis=0)
        self.y_pred_concat = pd.DataFrame(data=None, index=self.Y.index, columns=self.Y.columns)

        for fold in range(self.n_folds):
            self.y_pred_concat.loc[self.testing_genes[fold]] = self.y_pred[fold]
        t2 = time.time()
        print('Time elapsed: %.2f seconds'%(t2-t1))

    def grid_search(self, n_steps, params, fixed_params=dict(), verbose=False, stages=2, axis=0):
        """Perform grid search for hyperparameter optimization."""
        if verbose:
            print("stages remaining: " + str(stages))
        param_names = list(params.keys())
        n_params = len(param_names)
        
        # Create log-spaced parameter values
        param_steps = [np.log10(params[param_names[i]][1] / params[param_names[i]][0]) / (n_steps - 1) 
                      for i in range(n_params)]
        param_mins = [np.log10(params[param_names[i]][0]) for i in range(n_params)]
        
        # Track performance
        perf = list()
        best_perf = -1
        best_params = None
        
        # Progress bar
        if verbose:
            from tqdm import tqdm
            pbar = tqdm(range(n_steps ** n_params))
        else:
            pbar = range(n_steps ** n_params)
            
        # Evaluate all parameter combinations
        for i in pbar:
            step = [(i // (n_steps ** j)) % n_steps for j in range(n_params)]
            params_i = dict(
                [(param_names[j], 10 ** (param_mins[j] + step[j] * param_steps[j])) 
                 for j in range(n_params)])
            self.train(**params_i, **fixed_params)
            cor, gene_cor = self.evaluate()
            
            # Track best performing parameters
            if np.mean(cor) > best_perf:
                best_perf = np.mean(cor)
                best_params = params_i

        # Recursively refine search or return results
        if stages == 1 or best_params is None:
            if best_params is None:
                best_params = params_i
            train_params = {**best_params, **fixed_params}
            self.train(**train_params)
            return train_params
        else:
            # Narrow search range around best parameters
            params = dict([(param_names[j], 
                          [10 ** (np.log10(best_params[param_names[j]]) - param_steps[j]),
                           10 ** (np.log10(best_params[param_names[j]]) + param_steps[j])]) 
                         for j in range(n_params)])
            return self.grid_search(n_steps, params, fixed_params=fixed_params, 
                                  verbose=verbose, stages=stages-1, axis=axis)


class MultiRidge(ModelBase):
    """Multi-task Ridge regression model"""
    
    def __init__(self, args, **kwargs):
        super().__init__(args, **kwargs)

    def train(self, fold=0, lam=1):
        """Train ridge regression model for given fold."""
        if self.svdD[fold] is None:
            self.svdD[fold] = np.linalg.svd(self.D_train[fold], full_matrices=False)
            
        (ua, sa, va) = self.svdD[fold]
        # Compute ridge regression solution using SVD
        self.W[fold] = va.T.dot(
            np.diag(1/(lam*self.n_genes+sa**2)).dot(va)
        ).dot(self.D_train[fold].T.dot(self.Y_train[fold]))


class Stan(ModelBase):    
    def __init__(self, args, **kwargs):
        super().__init__(args, **kwargs)

    def train(self, fold=0, lam2=1, lam1=5):
        """Train STAN model for given fold."""
        if self.svdD[fold] is None:
            self.svdD[fold] = np.linalg.svd(self.D_train[fold], full_matrices=False)
        if self.svdK is None:
            # Compute SVD of spatial kernel matrix
            uk, sk, vk = np.linalg.svd(self.K, full_matrices=True)
            sk = np.concatenate((sk, [0]*(uk.shape[0]-len(sk))))
            self.svdK = (uk, sk, vk)

        [ud, sd, vd] = self.svdD[fold]
        [uk, sk, vk] = self.svdK
        D = self.D_train[fold]
        Y = self.Y_train[fold]

        # Compute spatially-regularized solution
        scale = np.divide(1, sd.reshape((-1,1))**2 + lam1*lam2/(lam1*sk**2+lam2).reshape(1,-1))
        W = vd.T@(np.multiply(scale, vd@D.T@Y@uk))@uk.T
        self.W[fold] = W