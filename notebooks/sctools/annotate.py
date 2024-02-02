import scanpy as sc
import pandas as pd


def majority_vote(adata, prediction_col, clustering_resolution, majority_col_name = None):
    # partly taken from celltypist
    key_added = f'leiden_scvi_{clustering_resolution}'
    sc.tl.leiden(
        adata, 
        key_added = key_added,
        resolution = clustering_resolution
    )
    clustering = adata.obs.pop(key_added)
    votes = pd.crosstab(adata.obs[prediction_col], clustering)
    majority = votes.idxmax(axis=0)
    majority = majority[clustering].reset_index()
    majority.index = adata.obs.index
    
    majority_col_name = majority_col_name if majority_col_name else 'majority_voting'
    colnames = ['clustering', majority_col_name]
    majority.columns = colnames
    majority[majority_col_name] = majority[majority_col_name].astype('category')
    
    for col in colnames:
        if col in adata.obs.columns:
            adata.obs.pop(col)
    
    adata.obs = adata.obs.join(majority)
    