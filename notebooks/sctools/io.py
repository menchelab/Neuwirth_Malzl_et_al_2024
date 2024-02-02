import gzip

import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np

from scipy.io import mmread, mmwrite


def read_sc_data(
    counts_file, 
    features_file,
    metadata_file
):
    data_dict = {}
    for key, filename in zip(
        ['counts', 'features', 'metadata'],
        [counts_file, features_file, metadata_file]
    ):
        if filename.endswith('gz'):
            open_file = lambda x: gzip.open(x, 'rt')
            
        else: 
            open_file = lambda x: open(x, 'r')
            
        with open_file(filename) as file:
            if key == 'counts':
                # transpose due to the way the data was exported to comply with Seurat
                # see also convert_to_raw.ipynb
                data = mmread(file).T.tocsr()
                
            elif key == 'metadata':
                data = pd.read_csv(
                    file,
                    sep = '\t',
                    index_col = 0
                )
            
            else:
                data = pd.DataFrame(
                    index = file.read().rstrip().split()
                )
            
        data_dict[key] = data

    adata = ad.AnnData(
        X = data_dict['counts'],
        obs = data_dict['metadata'],
        var = data_dict['features'],
        dtype = np.int64
    )
    return adata


def initialize_from_raw(file):
    adata = sc.read_h5ad(file)
    obs = adata.obs.copy()
    var = adata.raw.var.copy()
    adata = ad.AnnData(
        X = adata.raw.X,
        obs = obs,
        var = var,
        dtype = np.int64
    )
    return adata


def write_sc_data(adata, outfile_prefix, layer = None, obs_columns = None, var_columns = None):
    obs = adata.obs
    var = adata.var
    X = adata.layers[layer] if layer else adata.X
    if not isinstance(obs_columns, type(None)):
        obs = obs.loc[:, obs_columns]

    if not isinstance(var_columns, type(None)):
        var = var.loc[:, var_columns]

    mmwrite(
        '.'.join([outfile_prefix, layer if layer else 'X', 'mtx']),
        X
    )
    obs.to_csv(
        '.'.join([outfile_prefix, 'obs.tsv']),
        sep = '\t'
    )
    var.to_csv(
        '.'.join([outfile_prefix, 'var.tsv']),
        sep = '\t'
    )
