import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
import glob, os


def run_decoupler(adata, net):
    # Input: raw scRNA-seq data
    #        TF-gene dataframe
    adata.raw = adata
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    dc.run_ulm(mat=adata, net=net, source='source', target='target', weight='weight', verbose=True)
    dc.run_mlm(mat=adata, net=net, source='source', target='target', weight='weight', verbose=True)
    return adata


# from arboreto.algo import grnboost2
# from pyscenic.utils import modules_from_adjacencies, load_motifs
# from pyscenic.prune import prune2df, df2regulons
# from pyscenic.aucell import aucell


# def run_aucell(adata, net, dbs, motif_annotation_file):
#     # Input: raw count matrix
#     #        TF-gene dataframe (adjacencies)

#     if not os.path.exists('pyscenic'):
#         os.makedirs('pyscenic')

#     ex_matrix = adata.to_df()
#     genes = np.intersect1d(net['target'], ex_matrix.columns)
#     tfs = np.intersect1d(net['TF'], ex_matrix.columns)
#     ex_matrix = ex_matrix[list(set(genes).union(set(tfs)))]
#     net = net[net['TF'].isin(tfs)]
#     net = net[net['target'].isin(genes)]
#     adjacencies = net

#     modules = modules_from_adjacencies(adjacencies, ex_matrix, rho_mask_dropouts = True)
#     df = prune2df(dbs, modules, motif_annotation_file)

#     if not os.path.exists('pyscenic'):
#         os.makedirs('pyscenic')
#     df.to_csv('pyscenic/motifs.csv')
#     df_motifs = load_motifs('pyscenic/motifs.csv')

#     regulons = df2regulons(df_motifs)
#     auc_mtx = aucell(ex_matrix, regulons, noweights = True)

#     tfs = [x[:-3] for x in auc_mtx.columns]
#     auc_mtx.columns = tfs
#     return auc_mtx


# def run_pyscenic(adata, tf_names, dbs, motif_annotation_file,
#     run_grnboost2=True, adjacencies_file='pyscenic/adjacencies.csv'):
#     if not os.path.exists('pyscenic'):
#         os.makedirs('pyscenic')

#     ex_matrix = adata.to_df()
#     if run_grnboost2:
#         adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
#         adjacencies.to_csv(adjacencies_file, index=True)
#     else:
#         adjacencies = pd.read_csv(adjacencies_file, index_col=0)

#     modules = modules_from_adjacencies(adjacencies, ex_matrix, rho_mask_dropouts = True)
#     df = prune2df(dbs, modules, motif_annotation_file)

#     df.to_csv('pyscenic/motifs.csv')
#     df_motifs = load_motifs('pyscenic/motifs.csv')

#     regulons = df2regulons(df_motifs)
#     auc_mtx = aucell(ex_matrix, regulons, noweights = True)
#     tfs = [x[:-3] for x in auc_mtx.columns]
#     auc_mtx.columns = tfs
#     return auc_mtx