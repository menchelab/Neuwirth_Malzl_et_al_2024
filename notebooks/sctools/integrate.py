import scvi

import scanpy as sc

def filter_too_few_cell_batches(adata, batch_key, min_cells = 5):
    cells_per_batch = adata.obs.groupby(batch_key).count().iloc[:, 0]
    print(
        f'filtered the following batches < {min_cells} cells: ', 
        cells_per_batch[cells_per_batch < min_cells]
    )
    enough_cells = adata.obs[batch_key].apply(
        lambda x: cells_per_batch[x] >= min_cells
    )
    enough_cells_idx = enough_cells[enough_cells]
    return adata[enough_cells_idx].copy()


def integrate_data_scvi(
    adata, 
    batch_key, 
    categorical_covariate_keys = None,
    continuous_covariate_keys = None,
    use_highly_variable_genes = True,
    n_top_genes = 4000,
    use_gpu = True,
    max_epochs = None,
    train_size = 0.9,
    filter_small_batches = True
    
):
    if filter_small_batches:
        adata = filter_too_few_cell_batches(adata, batch_key)

    else:
        adata = adata.copy()

    adata.layers['counts'] = adata.X.copy()
    adata.raw = adata
    
    if use_highly_variable_genes:
        print('computing highly variable genes')
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes = n_top_genes,
            layer = 'counts',
            subset = True,
            flavor = 'seurat_v3',
        )
        
    scvi.model.SCVI.setup_anndata(
        adata,
        layer = 'counts',
        batch_key = batch_key,
        categorical_covariate_keys = categorical_covariate_keys,
        continuous_covariate_keys = continuous_covariate_keys
    )
    # non default parameters from scVI tutorial and scIB github
    # see https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html
    # and https://github.com/theislab/scib/blob/main/scib/integration.py
    model = scvi.model.SCVI(
        adata,
        n_layers = 2,
        n_latent = 30,
        gene_likelihood = 'nb'
    )
    model.train(
        use_gpu = use_gpu,
        max_epochs = max_epochs,
        train_size = train_size
    )
    adata.obsm['X_scvi'] = model.get_latent_representation()
    
    print('compute umap from scvi embedding')
    sc.pp.neighbors(
        adata,
        use_rep = 'X_scvi'
    )
    sc.tl.umap(
        adata
    )
    
    return {'data': adata, 'model': model}
