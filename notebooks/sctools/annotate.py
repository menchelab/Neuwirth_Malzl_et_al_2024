import scanpy as sc
import pandas as pd
import numpy as np


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


def annotate_adata_on_gene_hi_lo(adata, gene, resolution = 5):
    expression_values = adata[:, gene].X.toarray().flatten()
    threshold = np.median(expression_values)
    status_column = f'{gene.lower()}_status'
    adata.obs[status_column] = np.select(
        [expression_values < threshold, expression_values >= threshold],
        [f'{gene}_lo', f'{gene}_hi']
    )
    majority_vote(
        adata,
        status_column,
        resolution,
        f'{status_column}_majority_vote'
    )

def get_cells_per_nhood(nhoods_graph, row_index, nhood_center_cells):
    n_rows, n_cols = nhoods_graph.shape
    nhood_cells = {}
    for j, center_cell in enumerate(nhood_center_cells):
        nhood_neighbours = nhoods_graph[:, j]
        nhood_cells[center_cell] = row_index[
            (nhood_neighbours == 1).toarray().flatten()
        ]
    
    return nhood_cells


def assign_to_overrepresented_nhoods(adata, fdr = 0.1):
    cells = adata.obs.index
    nhood_centers = adata.obs.nhood_ixs_refined == 1
    nhood_data = adata.uns['nhood_adata'].obs.set_index('index_cell')
    cells_per_nhood = get_cells_per_nhood(
        adata.obsm['nhoods'],
        cells,
        cells[nhood_centers]
    )
    cells_in_overrepresented_nhoods = set()
    for center_cell, cells_in_nhood in cells_per_nhood.items():
        spatial_fdr = nhood_data.loc[center_cell, 'SpatialFDR']
        log2FC = nhood_data.loc[center_cell, 'logFC']
        if (spatial_fdr <= fdr) & (log2FC > 0):
            cells_in_overrepresented_nhoods.update(
                cells_in_nhood
            )
        
    adata.obs['nhood_annotation'] = 'else'
    adata.obs.loc[list(cells_in_overrepresented_nhoods), 'nhood_annotation'] = 'overrepresented'
