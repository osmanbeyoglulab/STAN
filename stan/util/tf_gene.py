import pandas as pd
import decoupler
import numpy as np
from anndata import AnnData
import scanpy as sc
from typing import Optional, Union, Literal, List
from pandas import DataFrame
import os.path
import os

def make_dir(dir):
    isExist = os.path.exists(dir)
    if not isExist:
       # Create a new directory because it does not exist
       os.makedirs(dir)

def tf_gene_cor(gex_df, tf_gene):
    genes=np.intersect1d(tf_gene['gene'].unique(), gex_df.columns)
    tfs=np.intersect1d(tf_gene['TF'].unique(), gex_df.columns)
    tf_gene=tf_gene.query('gene in @genes and TF in @tfs')

    genes=tf_gene['gene'].unique()
    tfs=tf_gene['TF'].unique()
    tf_gene=tf_gene.query('gene in @genes and TF in @tfs')

    print("n_genes: ", len(genes), "n_tfs: ", len(tfs))

    cors=dict()
    for tf in tfs:
        cors[tf]=gex_df[genes].corrwith(gex_df[tf])
    cors=pd.DataFrame(cors).reset_index().melt(id_vars='index')

    cors.columns=['gene','TF', 'cor']
    cors=cors[['TF', 'gene', 'cor']]

    tf_gene=tf_gene.pivot(index='gene', columns='TF', values='weight').reset_index().melt(id_vars='gene').fillna(0)[['TF', 'gene', 'value']]

    cors.sort_values(by=['TF', 'gene'])
    tf_gene.sort_values(by=['TF', 'gene'])

    cors['tf_gene']=tf_gene['value']
    cors["tf_gene"]=cors["tf_gene"].astype("category")

    pct=cors.pivot(index='gene', columns='TF', values='cor').rank(pct=True).reset_index().melt(id_vars='gene')[['TF', 'gene', 'value']]
    pct.sort_values(by=['TF', 'gene'])
    cors['pct']=pct['value']

    return cors


def tf_gene(
        data: AnnData,
        min_cells_proportion: float = 0.2,
        min_tfs_per_gene: int = 5,
        min_genes_per_tf: int = 10,
        max_proportion_genes_per_tf: float = 0.8,
        tf_gene_source: Union[Literal["dorothea", "hTFtarget"], DataFrame]="hTFtarget",
        tf_list: Union["humantfs", List, None]="humantfs",
        inplace: bool = True,
        source_dir: str = "data/tf_gene/"
        ) -> Union[AnnData, None]:

    '''
    Parameters
    -----------
    data
        An annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_cells_proportion: (between 0 and 1)
        Minimum proportion of cells expressed required for a gene to pass filtering
    min_tfs_per_gene:
        Minimum number of TFs that a target a gene required for that gene to pass filtering
    min_genes_per_tf:
        Minimum number of genes a TF must target to pass filtering
    max_proportion_genes_per_tf:
        Maximum percentage of genes a TF can target to pass filtering
    tf_gene_source:
        Source for tfs and target genes. Either the name of the database used or a pandas dataframe with columns 'TF' and 'gene'
    tf_list:
        List of TFs to limit analysis.  If "humantfs", then the list comes from TFs with known motifs in the HumanTFs database (http://humantfs.ccbr.utoronto.ca/allTFs.php)


    Returns
    -------
    Depending on `inplace`, returns a filtered AnnData object with less genes and TF-gene matrix in as adata.varm['tf_gene'], or None.

    '''

    adata=data if inplace else data.copy()

    tf_universe=adata.var_names.copy()

    sc.pp.filter_genes(adata,
                       min_cells=min_cells_proportion*adata.n_obs)


    tf_gene=make_tf_gene_matrix(source=tf_gene_source,
                                tf_list="humantfs",
                                gene_universe=adata.var_names,
                                tf_universe=tf_universe,
                                min_tfs_per_gene=min_tfs_per_gene,
                                min_genes_per_tf=min_genes_per_tf
                                )

    var_names=tf_gene.index
    tf_names=tf_gene.columns

    adata=adata[:, var_names]
    adata.varm['tf_gene']=tf_gene
    adata.uns['tf_names']=tf_names.to_list()

    if inplace==False:
        return adata




def make_tf_gene_matrix(source="hTFtarget",
                        tf_list="humantfs",
                        gene_universe=None,
                        tf_universe=None,
                        min_tfs_per_gene=10,
                        min_genes_per_tf=10,
                        tissue=None,
                        levels=['A', 'B', 'C'],
                        source_dir: str = "data/tf_gene/"
                        ):

    if type(source)== pd.DataFrame:
        tf_gene=source
    elif source=="dorothea":
        tf_gene=load_dorothea(levels=levels, source_dir=source_dir)
    elif source=="hTFtarget":
        tf_gene=load_hTFtarget(tissue=tissue, source_dir=source_dir)


    if tf_list=="humantfs":
        tf_list=load_humantfs(source_dir=source_dir)

    if tf_list is not None:
        tf_gene=tf_gene.query("TF in @tf_list")

    if tf_universe is not None:
        tf_gene=tf_gene.query("TF in @tf_universe") #just say >0

    if gene_universe is not None:
        tf_gene=tf_gene.query("gene in @gene_universe")

    tf_gene['values']=1
    tf_gene=tf_gene.pivot_table(index='gene', columns='TF', aggfunc='mean', values='values',fill_value=0)

    gene_sum=tf_gene.T.sum()
    tf_sum=tf_gene.sum()

    tf_gene.drop(tf_gene.index[gene_sum<min_tfs_per_gene], 0,inplace=True)
    tf_gene.drop(tf_gene.columns[tf_sum<min_genes_per_tf], 1,inplace=True)

    return tf_gene

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

def load_dorothea(levels=['A', 'B', 'C']):
    dorothea=decoupler.get_dorothea(levels=levels)
    dorothea=dorothea[["source", "target"]]
    dorothea.columns=["TF", "gene"]
    return dorothea

def load_humantfs(source_dir="data/tf_gene/"):
    fname=source_dir+"humantfs.csv"
    if os.path.isfile(fname):
        humantfs=pd.read_csv(fname)
    else:
        print('downloading humantfs database...')
        humantfs=pd.read_csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")
        make_dir(source_dir)
        humantfs.to_csv(fname)

    return humantfs.query("`Is TF?`=='Yes' and `TF assessment` == 'Known motif'")['HGNC symbol'].tolist()

def load_hTFtarget(tissue=None, source_dir="data/tf_gene/"):

    fname=source_dir+"hTFtaget.csv"
    if os.path.isfile(fname):
        htftarget=pd.read_csv(fname)
    else:
        print('downloading hTFtarget database...')
        hTFtarget_url="http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
        htftarget=pd.read_csv(hTFtarget_url,sep='\t')
        make_dir(source_dir)
        htftarget.to_csv(fname)

    if tissue is not None:
        include= [tissue in x for x in htftarget['tissue']]
        htftarget=htftarget.iloc[include]

    htftarget=htftarget.drop(columns="tissue")
    htftarget.columns=["TF", "gene"]

    return htftarget