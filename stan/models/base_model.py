from ..util.spatrafact_util import *
from ..pp import *
import numpy as np
from scipy.sparse.linalg import svds
from scipy.stats import pearsonr
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import FastICA
from sklearn.manifold import TSNE
from sklearn.preprocessing import QuantileTransformer

from copy import deepcopy
from tqdm.auto import tqdm
from copy import deepcopy
import pickle5 as pickle

class ModelBase:
    def __init__(self, adata, layer="dca", tf_gene=None, intercept=False):
        self.adata = adata
        self.intercept = intercept
        if tf_gene is None:
            tf_gene = adata.varm['tf_gene']

        self.tf_gene = tf_gene
        counts=adata.to_df(layer)

        self.Y = counts.T

        self.n_folds = np.max(adata.var['fold']) + 1

        self.genes = np.intersect1d(self.Y.index, self.tf_gene.index)
        self.n_genes = len(self.genes)

        self.tf_gene = self.tf_gene.loc[self.genes]
        self.Y = self.Y.loc[self.genes]

        self.n_spots = adata.n_obs

        self.training_genes = [adata.var.query("fold != @i").index for i in range(self.n_folds)]
        self.testing_genes = [adata.var.query("fold == @i").index for i in range(self.n_folds)]

        self.n_tfs = self.tf_gene.shape[1]

        self.D = self.tf_gene.to_numpy()

        self.D_train = [self.tf_gene.loc[self.training_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.D_test = [self.tf_gene.loc[self.testing_genes[i]].to_numpy() for i in range(self.n_folds)]

        self.Y_train = [self.Y.loc[self.training_genes[i]].to_numpy() for i in range(self.n_folds)]
        self.Y_test = [self.Y.loc[self.testing_genes[i]].to_numpy() for i in range(self.n_folds)]

        if intercept:
            self.D = np.concatenate((self.D, 1000 * np.ones((self.D.shape[0], 1))), axis=1)
            for i in range(self.n_folds):
                self.D_train[i] = np.concatenate((self.D_train[i], 1000 * np.ones((self.D_train[i].shape[0], 1))),
                                                 axis=1)
                self.D_test[i] = np.concatenate((self.D_test[i], 1000 * np.ones((self.D_test[i].shape[0], 1))), axis=1)

        self.kernel = adata.obsm.get("kernel")

        self.svdD=[]
        for i in range(self.n_folds):
            self.svdD.append(np.linalg.svd(self.D_train[i], full_matrices=False))

        self.svdK = None

        self.V = [None] * self.n_folds
        self.S = [None] * self.n_folds

        self.W = [None] * self.n_folds
        self.y_pred = [None] * self.n_folds
        self.tf_p_vals = None
        self.W_p_vals = None

        self.tf_names = self.tf_gene.columns.to_list()
        if self.intercept:
            self.tf_names.append('intercept')

    def results_to_adata(self, suffix="", outfile="sparafact.h5ad"):

        self.adata.obsm["tf_activity" + suffix] = pd.DataFrame(self.W_concat.T, columns=self.tf_names,
                                                               index=self.adata.obs_names)
        self.adata.layers["y_pred" + suffix] = self.y_pred_concat.T
        self.adata.layers["y_fit" + suffix] = (self.D.dot(self.W_concat)).T

        # # TODO: write this
        # if self.kernel is not None:
        #     uk, sk, vk = np.linalg.svd(self.kernel, full_matrices=False)
        #     proj = uk.dot(uk.T)
        #     pd.Series(np.sum((self.W_concat.dot(proj))**2, axis=1)/np.sum(self.W_concat**2, axis=1), index=colums)
        self.adata.uns["tf_spatial_score"] = self.tf_spatial_score()

        if outfile is not None:
            sc.write(outfile, self.adata)
    def get_tf(self):
        return pd.DataFrame(self.W_concat.T, columns=self.tf_names,index=self.adata.obs_names)
    def evaluate(self, fold=0, return_string=False, gene_set="testing", mask=None):
        if fold == -1:
            y_pred = np.asarray(self.y_pred_concat)
            Y = np.asarray(self.Y)
        else:
            W = self.W[fold]
            if gene_set == "testing":
                Y = self.Y_test[fold]
                D = self.D_test[fold]
            elif gene_set == "training":
                Y = self.Y_train[fold]
                D = self.D_train[fold]
            y_pred = D.dot(W)

        if mask is None:
            cor = [pearsonr(y_pred[:, spot], Y[:, spot])[0] for spot in range(self.n_spots)]
            gene_cor = [pearsonr(y_pred[i, :].T, Y[i, :].T)[0] for i in range(Y.shape[0])]

        if return_string:
            return "Sample Cor:" + str(round(np.nanmedian(cor), 3)) + " Gene Cor: " + str(round(np.nanmedian(gene_cor), 3))
        else:
            return cor, gene_cor

    def fit(self, grid_search_params=None, fixed_params=dict(), verbose=False, stages=5, n_steps=5, axis=0):
        if grid_search_params is None:
            self.params = fixed_params
        else:
            self.params = self.grid_search(n_steps, grid_search_params, fixed_params=fixed_params, verbose=verbose, stages=stages,
                                      axis=axis)
        params=self.params
        for fold in range(self.n_folds):
            self.train(fold=fold, **params)
            self.y_pred[fold] = self.D_test[fold].dot(self.W[fold])
            if verbose:
                print("Fold " + str(fold) + " Fit Accuracy        " + self.evaluate(fold=fold, return_string=True,
                                                                                    gene_set="training"))
                print("Fold " + str(fold) + " Prediction Accuracy " + self.evaluate(fold=fold, return_string=True))

        self.W_concat = np.mean(self.W, axis=0)
        self.y_pred_concat = pd.DataFrame(data=None, index=self.Y.index, columns=self.Y.columns)

        for fold in range(self.n_folds):
            self.y_pred_concat.loc[self.testing_genes[fold]] = self.y_pred[fold]

    def grid_search(self, n_steps, params, fixed_params=dict(), verbose=False, stages=2, axis=0):
        if verbose:
            print("stages remaining: " + str(stages))
        param_names = list(params.keys())
        n_params = len(param_names)
        param_steps = [np.log10(params[param_names[i]][1] / params[param_names[i]][0]) / (n_steps - 1) for i in
                       range(n_params)]
        param_mins = [np.log10(params[param_names[i]][0]) for i in range(n_params)]
        perf = list()
        best_perf = -1
        best_params= None
        if verbose:
            pbar = tqdm(range(n_steps ** n_params))
        else:
            pbar = range(n_steps ** n_params)
        for i in pbar:
            step = [(i // (n_steps ** j)) % n_steps for j in range(n_params)]
            params_i = dict(
                [(param_names[j], 10 ** (param_mins[j] + step[j] * param_steps[j])) for j in range(n_params)])
            self.train(**params_i, **fixed_params)
            if np.mean(self.evaluate()[axis]) > best_perf:
                best_perf = np.mean(self.evaluate()[axis])
                best_params = params_i

        if stages == 1 or best_params is None:
            if best_params is None:
                best_params = params_i
            train_params = {**best_params, **fixed_params}
            self.train(**train_params)
            return train_params
        else:
            params = dict([(param_names[j], [10 ** (np.log10(best_params[param_names[j]]) - param_steps[j]),
                                             10 ** (np.log10(best_params[param_names[j]]) + param_steps[j])]) for j in
                           range(n_params)])
            return self.grid_search(n_steps, params, fixed_params=fixed_params, verbose=verbose, stages=stages - 1,
                                    axis=axis)

    def permutation_test(self, n_permutations=100):

        perm_copy = deepcopy(self)

        W_p_vals_np = np.zeros(self.W_concat.shape)
        tf_p_vals_np = np.zeros(self.W_concat.shape[0])
        # self.tf_score_p_val=np.zeros(tf_scores.shape)

        for i in tqdm(range(n_permutations)):

            Y_perm = perm_copy.Y.to_numpy()
            Y_perm = Y_perm[np.random.permutation(perm_copy.Y.shape[0]), :]
            # Y_perm=Y_perm[:,np.random.permutation(perm_copy.Y.shape[1])]

            perm_copy.Y = pd.DataFrame(Y_perm, index=perm_copy.Y.index, columns=perm_copy.Y.columns)

            perm_copy.Y_train = [perm_copy.Y.loc[perm_copy.training_genes[i]].to_numpy() for i in
                                 range(perm_copy.n_folds)]
            perm_copy.Y_test = [perm_copy.Y.loc[perm_copy.testing_genes[i]].to_numpy() for i in
                                range(perm_copy.n_folds)]

            for fold in range(perm_copy.n_folds):
                perm_copy.train(fold=fold, **perm_copy.params)

            W_perm = np.mean(perm_copy.W, axis=0)

            W_p_vals_np += (1 / n_permutations) * np.less(np.abs(self.W_concat), np.abs(W_perm))
            tf_p_vals_np += (1 / n_permutations) * np.less(np.linalg.norm(self.W_concat, axis=1),
                                                           np.linalg.norm(W_perm, axis=1))

        self.W_p_vals = pd.DataFrame(W_p_vals_np, columns=self.adata.obs_names, index=self.tf_names).T
        self.tf_p_vals = pd.Series(tf_p_vals_np, index=self.tf_names)

        return W_p_vals_np
