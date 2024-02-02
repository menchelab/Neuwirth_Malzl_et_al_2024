import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
import numpy as np


mpl.rcParams['pdf.fonttype'] = 42


def plot_leiden_clusterings(
    data_dict, 
    resolutions, 
    data_key = None, 
    legend_loc = 'on_data', 
    panelheight = 5, 
    panelwidth = 6, 
    subsample = 1, 
    size = None
):
    fig, axs = plt.subplots(
        len(data_dict), 
        len(resolutions)
    )
    
    for i,(k, d) in enumerate(data_dict.items()):
        if not data_key:
            tmp = d.copy()
        
        else:
            tmp = d[data_key].copy()
        
        idx = np.random.choice(
            tmp.obs.index,
            size = int(tmp.shape[0] * subsample),
            replace = False
        )
        
        one_dimensional = len(data_dict) == 1
        for ax, resolution in zip(axs if one_dimensional else axs[i, :], resolutions):
            sc.tl.leiden(
                tmp, 
                key_added = f'leiden_scvi_{resolution}',
                resolution = resolution
            )
            sc.pl.umap(
                tmp[idx],
                color = f'leiden_scvi_{resolution}',
                frameon = False,
                show = False,
                ax = ax,
                size = size,
                legend_loc = legend_loc
            )

        del tmp
    
    fig.set_figwidth(panelwidth * len(resolutions))
    fig.set_figheight(panelheight * len(data_dict))
    fig.tight_layout()
    
    return fig, axs


def plot_integration_results(data_dict, color_keys, params_list = None, data_key = None, panelheight = 5, panelwidth = 6, subsample = 1):
    fig, axs = plt.subplots(len(color_keys), len(data_dict))
    for i, (k, d) in enumerate(data_dict.items()):
        if data_key:
            data = d[data_key]
            
        else:
            data = d
            
        if not params_list:
            params_list = [{} for i in range(len(color_keys))]
            
        idx = np.random.choice(
            data.obs.index,
            size = int(data.shape[0] * subsample),
            replace = False
        )
        one_dimensional = len(data_dict) == 1
        for ax, color_key, kwargs in zip(axs if one_dimensional else axs[:, i], color_keys, params_list):
            sc.pl.umap(
                data[idx],
                color = color_key,
                show = False,
                frameon = False,
                ax = ax,
                **kwargs
            )
        
        top_ax = axs[0] if one_dimensional else axs[0, i]
        top_ax.set_title(k)
        
    fig.set_figwidth(panelwidth * len(data_dict))
    fig.set_figheight(panelheight * len(color_keys))
    fig.tight_layout()
    
    return fig, axs


def plot_clustering_and_expression(
    data_dict, 
    cluster_keys,
    expression_keys, 
    params_list = None, 
    data_key = None,     
    panelheight = 5, 
    panelwidth = 6, 
    subsample = 1
):
    fig, axs = plt.subplots(len(expression_keys) + 1, len(data_dict))
    for i, (k, d) in enumerate(data_dict.items()):
        color_keys = [cluster_keys[k]] + expression_keys
        if data_key:
            data = d[data_key]
        
        else:
            data = d
        
        if not params_list:
            params_list = [{} for i in range(len(color_keys))]
            
        idx = np.random.choice(
            data.obs.index,
            size = int(data.shape[0] * subsample),
            replace = False
        )
            
        one_dimensional = len(data_dict) == 1
        for ax, color_key, kwargs in zip(axs if one_dimensional else axs[:, i], color_keys, params_list):
            sc.pl.umap(
                data[idx],
                color = color_key,
                frameon = False,
                show = False,
                ax = ax,
                **kwargs
            )
        
        top_ax = axs[0, i] if len(axs.shape) > 1 else axs[0]
        top_ax.set_title(k)
    
    fig.set_figwidth(panelwidth * len(data_dict))
    fig.set_figheight(panelheight * len(color_keys))
    fig.tight_layout()
    
    return fig, axs


def raw_data_umap(adata, color, nhvg = 4000, **kwargs):
    tmp = adata.copy()
    sc.pp.normalize_total(
        tmp, 
        target_sum = 1e4
    )
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(
        tmp,
        n_top_genes = 4000,
        layer = "counts",
        flavor = "seurat_v3",
    )
    sc.pp.pca(
        tmp, 
        n_comps = 40, 
        svd_solver = 'arpack',
        use_highly_variable = True
    )
    sc.pp.neighbors(
        tmp,
        use_rep = 'X_pca'
    )
    sc.tl.umap(tmp)

    sc.pl.umap(
        tmp,
        color = color,
        frameon = False,
        **kwargs
    )
    del tmp
