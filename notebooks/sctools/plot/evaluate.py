import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
import seaborn as sns
import numpy as np


mpl.rcParams['pdf.fonttype'] = 42


def plot_expression_histogram(adata, layer = 'counts', bins = 50):
    adata = adata.copy()

    if layer:
        adata.X = adata.layers['counts'].copy()
    
    sc.pp.normalize_total(
        adata,
        target_sum = 1e4
    )
    sc.pp.log1p(adata)

    fig, ax = plt.subplots()
    x = adata[:, 'SAT1'].X.toarray().flatten()
    thres = np.median(x)
    hue = ['SAT1_hi' if xi > thres else 'SAT1_lo' for xi in x]
    sns.histplot(
        x = x,
        bins = bins,
        hue = hue,
        ax = ax,
        multiple = 'stack',
        palette = 'Set2'
    )
    ax.set_title('SAT1 expression')
    ax.set_xlabel('log(cpm)')
    ax.axvline(thres, ls = '--', c = 'grey')

    fig.set_figheight(5)
    fig.set_figwidth(10)
    fig.tight_layout()
    return fig


def annotate_enrichment(x, spatial_fdr_threshold = 0.25):
    if x['SpatialFDR'] < spatial_fdr_threshold:
        return 'enriched' if x['logFC'] > 0 else 'depleted'
    
    else:
        return 'not_significant'


def plot_nhood_violin(adata, spatial_fdr_threshold = 0.25, ax = None):
    if isinstance(ax, type(None)):
        fig, ax = plt.subplots()
        set_figure_extents = True
    
    else:
        set_figure_extents = False

    df = adata.uns['nhood_adata'].obs.copy()
    #df['nhood_size'] = np.array(adata.uns['nhood_adata'].X.sum(1)).flatten()
    df['enriched'] = df[['logFC', 'SpatialFDR']].apply(
        annotate_enrichment, 
        axis = 1,
        spatial_fdr_threshold = spatial_fdr_threshold
    )
    sns.violinplot(
        y = 'nhood_annotation',
        x = 'logFC',
        data = df,
        color = '#f2f2f2'
    )
    sns.stripplot(
        y = 'nhood_annotation',
        x = 'logFC',
        data = df,
        hue = 'enriched',
        palette = {
            'enriched': '#f44336',
            'not_significant': '#bcbcbc',
            'depleted': '#6fa8dc'
        },
        ax = ax,
        edgecolor = 'k',
        linewidth = 0.5
    )
    ax.set_ylabel('SAT1 status', fontsize = 20)
    ax.set_xlabel('logFC', fontsize = 20)
    ax.set_yticklabels(
        ['SAT1 low', 'SAT1 high']
    )
    ax.tick_params(
        labelsize = 15
    )
    ax.legend(
        loc='upper left', 
        bbox_to_anchor=(1, 1), 
        frameon=False,
        fontsize = 15
    )

    if set_figure_extents:
        fig.set_figwidth(10)
        fig.set_figheight(5)
        fig.tight_layout()

    return ax
