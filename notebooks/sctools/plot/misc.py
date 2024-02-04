import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
import seaborn as sns


mpl.rcParams['pdf.fonttype'] = 42


def generate_and_save_for_figure(
    adata,
    columns_to_plot,
    scatter_save_file,
    legend_save_file,
    scatter_dpi = 1000,
    scatter_size = (5, 5),
    size = 10,
    edgecolor = 'none',
    linewidths = 1.5
):
    for col, (palette, vmax) in columns_to_plot.items():
        # generate figure and remove legend
        fig, ax = plt.subplots()
        width, height = scatter_size
        colorbar = False
        
        # if column is not in obs it is a gene
        if col not in adata.obs.columns:
            color_map = palette
            palette = None
            colorbar = True
        
        elif palette is None:
            palette = None
            color_map = None
        
        elif isinstance(palette, str):
            ncolors = adata.obs[col].nunique()
            palette = sns.color_palette(palette, ncolors)
            color_map = None
        
        elif isinstance(palette, dict):
            color_map = None
        
        sc.pl.umap(
            adata,
            color = col,
            frameon = False,
            palette = palette,
            color_map = color_map,
            colorbar_loc = None, # no colorbar to ensure w/h ratio
            size = size,
            vmax = vmax,
            ax = ax,
            show = False,
            edgecolor = edgecolor,
            linewidths = linewidths
        )
            
        handles, labels = ax.get_legend_handles_labels()
        ax.legend().remove()

        ax.title.set_visible(False)
        
        fig.set_figwidth(width)
        fig.set_figheight(height)
        fig.tight_layout()
        plt.show()
        fig.savefig(scatter_save_file.format(col), dpi = scatter_dpi)
            
        
        if colorbar:
            plt.close(fig)
            fig, ax = plt.subplots()
            sc.pl.umap(
                adata,
                color = col,
                frameon = False,
                palette = palette,
                color_map = color_map,
                colorbar_loc = 'right',
                size = size,
                vmax = vmax,
                ax = ax,
                show = False
            )
        
            # remove scatters
            ax.collections.clear()
            fig.savefig(legend_save_file.format(col))
            plt.close(fig)
        
        if not colorbar:
            # remove scatters
            ax.collections.clear()
            ax.legend(handles, labels)
            fig.savefig(legend_save_file.format(col))
            plt.close(fig)
            