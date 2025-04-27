import anndata as ann
import json, os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from glob import glob
from matplotlib.image import imread
from pathlib import Path


def read_visium_sge(sample_id="V1_Human_Lymph_Node", min_cells=5, min_counts=5000):
    fpath = f'data/{sample_id}'
    if os.path.exists(f'{fpath}/filtered_feature_bc_matrix.h5'):
        adata = sq.read.visium(f'data/{sample_id}')
    else:
        adata = sc.datasets.visium_sge(sample_id)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata.layers['raw'] = adata.X.copy()
    return adata


def prepare_breast_wu(path, sample):
    matrix_path = f'{path}/filtered_count_matrices/'
    for f in glob(matrix_path + sample + "_filtered_count_matrix/*.gz"):
        os.rename(f, f[0:-3])
    adata = sc.read_mtx(matrix_path + sample + "_filtered_count_matrix/matrix.mtx").T
    with open(matrix_path + sample + "_filtered_count_matrix/barcodes.tsv", "r") as f:
        adata.obs_names = f.read().split("\n")[0:-1]
    with open(matrix_path + sample + "_filtered_count_matrix/features.tsv", "r") as f:
        adata.var_names = f.read().split("\n")[0:-1]

    adata.uns["spatial"] = dict()
    adata.uns["spatial"][sample] = dict()
    
    spatial_path = '{}/spatial/{}_spatial/'.format(path, sample)
    files = dict(
        tissue_positions_file = spatial_path + 'tissue_positions_list.csv',
        scalefactors_json_file = spatial_path + 'scalefactors_json.json',
        hires_image = spatial_path + 'tissue_hires_image.png',
        lowres_image = spatial_path + 'tissue_lowres_image.png',)

    adata.uns["spatial"][sample]['images'] = dict()
    for res in ['hires', 'lowres']:
        adata.uns["spatial"][sample]['images'][res] = imread(str(files[f'{res}_image']))

    # read json scalefactors
    adata.uns["spatial"][sample]['scalefactors'] = json.loads(Path(files['scalefactors_json_file']).read_bytes())

    # read coordinates
    positions = pd.read_csv(files['tissue_positions_file'], header=None)
    positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
    positions.index = positions['barcode']

    adata.obs = adata.obs.join(positions, how="left")
    adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()
    adata.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)
    
    metadata = pd.read_csv("{}/metadata/{}_metadata.csv".format(path, sample), index_col=0)
    adata.obs["subtype"] = metadata["subtype"]
    adata.obs["pathology"] = metadata["Classification"]
    adata.obs["sample"] = sample
    adata.obs["replicate"] = sample
    adata.obs["ER"] = [x=="ER" for x in adata.obs["subtype"]]
    adata.obs["HER2"] = False # all samples are HER2 - in this dataset.
    adata.obs["PR"] = adata.obs["ER"] # All ER+ samples are also PR+ in this dataset.
    adata.layers["raw_counts"]=adata.X
    adata.obs_names = sample+"_"+adata.obs_names
    return adata


def read_breast_wu(path2h5ad, min_cells=5, min_counts=500):
    # path2h5ad = path to h5ad file
    adata = sc.read_h5ad(path2h5ad)
    adata = adata[adata.obs['pathology'] != "Artefact"]
    adata = adata[adata.obs['pathology'] != "Uncertain"]
    adata = adata[adata.obs['pathology'] != np.nan]
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata.layers['raw'] = adata.X
    return adata


def read_cytassist(path, min_cells=5, min_counts=500, remove_isotype=True):
    data = sc.read_10x_h5(f"{path}/filtered_feature_bc_matrix.h5", gex_only=False)
    visium_ = sc.read_visium(path=path) 

    data.uns['spatial'] = visium_.uns['spatial'] 
    data.obsm['spatial'] = visium_.obsm['spatial']
    data.obsm['spatial'] = data.obsm['spatial'].astype(float)
    data.var_names_make_unique()
  
    adata = data[:, data.var.feature_types=='Gene Expression']
    pdata = data[:, data.var.feature_types=='Antibody Capture']

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
  
    if remove_isotype:
        pdata.var["isotype_control"] = (pdata.var_names.str.startswith("mouse_") \
        | pdata.var_names.str.startswith("rat_") \
        | pdata.var_names.str.startswith("HLA_")\
        | pdata.var_names.str.startswith("mouse.") \
        | pdata.var_names.str.startswith("rat.") \
        | pdata.var_names.str.startswith("HLA."))
 
    pdata = pdata[:, pdata.var.isotype_control==False]
    pdata = pdata[adata.obs_names,:]
    adata.layers['raw'] = adata.X.copy()
    pdata.layers['raw'] = pdata.X.copy()
    return adata, pdata

# read the dataset provided with 10x-Genomics-formatted mtx directory and spatial info. The directory should be named as: "filtered_feature_bc_matrix"

def read_cytassist_mtx(path, min_cells=5, min_counts=500, remove_isotype=True):  
    adata = sc.read_10x_mtx(f"{path}/filtered_feature_bc_matrix")
    adata.uns["spatial"] = dict()
    adata.uns["spatial"]["Visium"] = dict()
    
    spatial_path = f"{path}/spatial/"
    files = dict(
        tissue_positions_file = spatial_path + 'tissue_positions_list.csv',
        scalefactors_json_file = spatial_path + 'scalefactors_json.json',
        hires_image = spatial_path + 'tissue_hires_image.png',
        lowres_image = spatial_path + 'tissue_lowres_image.png',)

    adata.uns["spatial"]["Visium"]['images'] = dict()
    for res in ['hires', 'lowres']:
        if os.path.isfile(files[f'{res}_image']):
            adata.uns["spatial"]["Visium"]['images'][res] = imread(str(files[f'{res}_image']))

    # read json scalefactors
    adata.uns["spatial"]["Visium"]['scalefactors'] = json.loads(Path(files['scalefactors_json_file']).read_bytes())

    # read coordinates
    positions = pd.read_csv(files['tissue_positions_file'], header=None, index_col=0)
    positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']

    adata.obsm['spatial'] = positions.loc[adata.obs.index, ['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata.layers['raw'] = adata.X.copy()
    return adata

# read the dataset provided with raw counts dataframe and spatial info. The raw counts file should be named as: "raw_counts.csv"
def read_cytassist_csv(path, min_cells=5, min_counts=500, remove_isotype=True):
    df = pd.read_csv(f"{path}/raw_counts.csv", index_col=0)
    # the dataframe is gene X cells, and cell names might be modified by R code from "-" character to ".". 
    # Here we need to reverse back to match original name showed in position file
    df.columns = [ item.replace(".","-") for item in df.columns ]
    adata = ann.AnnData(df.T)
    
    adata.uns["spatial"] = dict()
    adata.uns["spatial"]["Visium"] = dict()
    
    spatial_path = f"{path}/spatial/"
    files = dict(
        tissue_positions_file = spatial_path + 'tissue_positions_list.csv',
        scalefactors_json_file = spatial_path + 'scalefactors_json.json',
        hires_image = spatial_path + 'tissue_hires_image.png',
        lowres_image = spatial_path + 'tissue_lowres_image.png',)

    adata.uns["spatial"]["Visium"]['images'] = dict()
    for res in ['hires', 'lowres']:
        if os.path.isfile(files[f'{res}_image']):
            adata.uns["spatial"]["Visium"]['images'][res] = imread(str(files[f'{res}_image']))

    # read json scalefactors
    adata.uns["spatial"]["Visium"]['scalefactors'] = json.loads(Path(files['scalefactors_json_file']).read_bytes())

    # read coordinates
    positions = pd.read_csv(files['tissue_positions_file'], header=None, index_col=0)
    positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']

    adata.obsm['spatial'] = positions.loc[adata.obs.index, ['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata.layers['raw'] = adata.X.copy()
    return adata


