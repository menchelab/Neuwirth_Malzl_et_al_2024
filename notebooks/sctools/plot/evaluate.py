import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd

from matplotlib import gridspec
from matplotlib.patches import Rectangle


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


def subplot_from_gridspec(
    fig, 
    gs, 
    row, 
    col, 
    show_xticks = False, 
    show_yticks = False, 
    show_spines = True
):
    ax = fig.add_subplot(gs[row, col])
    
    if not show_xticks:
        ax.set_xticks([])
        
    if not show_yticks:
        ax.set_yticks([])
        
    if not show_spines:
        for pos in ['top', 'bottom', 'left', 'right']:
            ax.spines[pos].set_visible(False)
        
    return ax


def setup_figure(width, height, x_groupby = [], y_groupby = []):
    fig = plt.figure(figsize = (width, height))
    heatmap_width = 8
    annotation_width = 0.2
    colorbar_width = 0.2
    legend_width = 0.2 * (len(x_groupby) + len(y_groupby)) * 2
    width_ratios = [annotation_width] * len(y_groupby) + [heatmap_width, colorbar_width]
    
    annotations = any(groupby for groupby in [x_groupby, y_groupby])
    if not annotations:
        height_ratios = [heatmap_width]
        heatmap_ax_row = 0
        
    else:
        height_ratios = [legend_width, heatmap_width] + [annotation_width] * len(x_groupby)
        heatmap_ax_row = 1
    
    gs = gridspec.GridSpec(
        nrows = len(height_ratios),
        ncols = len(width_ratios),
        width_ratios = width_ratios,
        height_ratios = height_ratios,
        wspace = 0.15 / width,
        hspace = 0.15 / height,
    )
    
    xticks_on_heatmap, yticks_on_heatmap = False, False
    if len(height_ratios) == 1:
        xticks_on_heatmap = True
        
    if len(width_ratios) == 2:
        yticks_on_heatmap = True
    
    axs = {
        'heatmap': subplot_from_gridspec(fig, gs, heatmap_ax_row, -2, xticks_on_heatmap, yticks_on_heatmap, False),
        'colorbar': subplot_from_gridspec(fig, gs, heatmap_ax_row, -1, False, False, False)
    }
    
    if annotations:
        axs['legend'] = subplot_from_gridspec(fig, gs, 0, -2, False, False, False)
    
    for j, group in enumerate(y_groupby):
        show_yticks = False
        if j == 0:
            show_yticks = True
            
        axs[group] = subplot_from_gridspec(fig, gs, heatmap_ax_row, j, False, show_yticks, False)
        
        if show_yticks:
            axs['y_axis'] = axs[group]
        
    for i, group in enumerate(x_groupby, heatmap_ax_row + 1):
        show_xticks = False
        if i == len(x_groupby) + heatmap_ax_row:
            show_xticks = True
            
        axs[group] = subplot_from_gridspec(fig, gs, i, -2, show_xticks, False, False)
        
        if show_xticks:
            axs['x_axis'] = axs[group]
        
    return fig, axs


def group_data(data, groupby):
    grouped_data = pd.concat(
        # avoid warning for single grouper by unpacking
        [g for _, g in data.groupby(groupby[0] if len(groupby) == 1 else groupby)]
    )
    group_info = {
        grouper: grouped_data.pop(grouper) for grouper in groupby
    }
    return grouped_data, group_info


def group_heatmap_data(data, x_groupby, y_groupby, x_group_df, y_group_df):
    data = data.merge(
        x_group_df,
        how = 'left',
        left_index = True,
        right_index = True
    )
    x_grouped_data, x_group_info = group_data(data, x_groupby)

    x_grouped_data = x_grouped_data.T.merge(
        y_group_df,
        how = 'left',
        left_index = True,
        right_index = True
    )
    xy_grouped_data, y_group_info = group_data(x_grouped_data, y_groupby)
    
    return xy_grouped_data, x_group_info, y_group_info


def groups_to_colors(group_info, palette_name):
    groups = sorted(group_info.unique())
    palette = sns.color_palette(palette_name, len(groups))
    colormap = {
        group: palette[i] for i, group in enumerate(groups)
    }
    return colormap


def group_info_to_colors(group_info, y_grouper, palette):
    if isinstance(palette, str):
        colormap = groups_to_colors(group_info, palette)
    
    else:
        colormap = palette

    color_array = np.zeros(shape = (1, len(group_info), 3))
    for i, (_, group_label) in enumerate(group_info.items()):
        color_array[0, i, :] = colormap[group_label]
    
    return np.swapaxes(color_array, 0, 1) if y_grouper else color_array


def add_group_annotation(
    group_info_dict, 
    axs, 
    group_palettes, 
    y_grouper = False, 
    plot_first_group_borders = True
):
    heatmap = axs['heatmap']
    for i, (group, group_info) in enumerate(group_info_dict.items()):
        ax = axs[group]
        # using pcolormesh here due to weird change of color
        # when importing in illustrator
        color_array = group_info_to_colors(
            group_info, 
            y_grouper, 
            palette = group_palettes[group]
        )
        ax.pcolormesh(color_array[::-1])
        
        if i == 0 and plot_first_group_borders:
            border = group_info.value_counts(sort = False).iloc[0] - 0.5
            line_drawer = heatmap.axhline if y_grouper else heatmap.axvline
            line_drawer(border, c = 'white', lw = 1)
        
        
def add_xgroup_annotation(group_info_dict, axs, x_group_palettes):
    add_group_annotation(group_info_dict, axs, x_group_palettes)
    

def add_ygroup_annotation(group_info_dict, axs, y_group_palettes):
    add_group_annotation(group_info_dict, axs, y_group_palettes, True)
        

def swap(a, b):
    tmp = a
    a = b
    b = tmp
    return a, b


def get_text_extents(text, renderer, inverse_ax_transform):
    bb = text.get_window_extent(renderer = renderer)
    transformed_bb = inverse_ax_transform.transform_bbox(bb)
    width, height = transformed_bb.width, transformed_bb.height
    x, y = transformed_bb.x0, transformed_bb.y0
    return x, y, width, height


def get_rect_extents(rect, inverse_ax_transform):
    transformed_bb = inverse_ax_transform.transform_bbox(
        rect.get_extents()
    )
    width, height = transformed_bb.width, transformed_bb.height
    x, y = transformed_bb.x0, transformed_bb.y0
    return x, y, width, height


def add_legend(
    fig, 
    ax, 
    x_group_info, 
    y_group_info, 
    x_group_palettes, 
    y_group_palettes, 
    patch_dim_ratio = 1
):   
    renderer = fig.canvas.get_renderer()
    inverse_ax_transform = ax.transData.inverted()
    vertical_spacing = 0.01
    legend_row_space = 0.1
    horizontal_spacing = 0.01
    patch_dim = 0.03
    legend_ypos = 0
    x_extent = 0
    for axis_label, group_info, group_palettes in zip(
        ['x-groups', 'y-groups'], 
        [x_group_info, y_group_info],
        [x_group_palettes, y_group_palettes]
    ):
        group_colors = {}
        for g, g_info in group_info.items():
            palette = group_palettes[g]
            if isinstance(palette, str):
                group_colors[g] = groups_to_colors(g_info, palette)

            else:
                group_colors[g] = palette

        n_legend_rows = len(group_info) * 2 # * 2 because group label is on top
        group_legend_height = legend_row_space * n_legend_rows + vertical_spacing * (n_legend_rows - 1) + legend_ypos
        alignment_ypos = group_legend_height - 1.5 * legend_row_space - vertical_spacing # 1.5 to put it in second row
        axis_label_text = ax.text(
            0.01, 
            alignment_ypos, 
            axis_label + ':', 
            fontsize = 8,
            va = 'center'
        )
        axis_label_x, _, axis_label_width, _ = get_text_extents(
            axis_label_text, 
            renderer,
            inverse_ax_transform
        )
        
        start_xpos = axis_label_x + axis_label_width + horizontal_spacing * 2
        current_group_label_ypos = group_legend_height - legend_row_space / 2
        # this actually needs to come from group_info
        # to make sure that we don't show groups that are not
        # in the group_info when using dict palettes
        for group_label, g_colors in group_colors.items():
            _ = ax.text(
                start_xpos,
                current_group_label_ypos,
                group_label,
                fontsize = 6,
                va = 'center'
            )
            current_xpos = start_xpos
            for label, color in g_colors.items():
                rectangle = Rectangle(
                    (current_xpos, alignment_ypos - patch_dim / 2), # va = 'center'
                    patch_dim, 
                    patch_dim / patch_dim_ratio,
                    color = color
                )
                rectangle = ax.add_patch(rectangle)
                _, _, rectangle_width, _ = get_rect_extents(
                    rectangle,
                    inverse_ax_transform
                )
                current_xpos = current_xpos + rectangle_width + horizontal_spacing
                
                label_text = ax.text(
                    current_xpos,
                    alignment_ypos,
                    label,
                    fontsize = 6,
                    va = 'center'
                )
                _, _, label_width, _= get_text_extents(
                    label_text,
                    renderer,
                    inverse_ax_transform
                )
                current_xpos = current_xpos + label_width + horizontal_spacing
            
            if current_xpos > x_extent:
                x_extent = current_xpos
            
            current_group_label_ypos = current_group_label_ypos - 2 * (legend_row_space + vertical_spacing)
            alignment_ypos = alignment_ypos - 2 * (legend_row_space + vertical_spacing)
                
        legend_ypos = legend_ypos + group_legend_height + vertical_spacing
        
    return group_legend_height, x_extent
    

def set_axes_extents(
    fig, 
    ax, 
    x_group_info, 
    y_group_info, 
    x_group_palettes, 
    y_group_palettes, 
    patch_dim_ratio = 1
):
    ymax, xmax = add_legend(
        fig,
        ax,
        x_group_info,
        y_group_info,
        x_group_palettes,
        y_group_palettes,
        patch_dim_ratio
    )
    
    ax.clear()
    
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_yticks([])
    ax.set_xticks([])
        
        
def grouped_heatmap(
    adata, 
    x_groupby, 
    y_groupby, 
    x_group_palettes,
    y_group_palettes,
    cmap, 
    vmin = None, 
    vmax = None, 
    show_var_labels = False, 
    swap_axes = False,
    figwidth = 10,
    figheight = 10
):
    data = adata.to_df()
    obs = adata.obs.loc[:, x_groupby]
    var = adata.var.loc[:, y_groupby]
    xy_grouped_data, x_group_info, y_group_info = group_heatmap_data(
        data, 
        x_groupby, 
        y_groupby, 
        obs, 
        var
    )
    
    if swap_axes:
        xy_grouped_data = xy_grouped_data.T
        x_group_info, y_group_info = swap(x_group_info, y_group_info)
        x_groupby, y_groupby = swap(x_groupby, y_groupby)
        x_group_palettes, y_group_palettes = swap(x_group_palettes, y_group_palettes)
        width_min = xy_grouped_data.shape[1] * 0.2 
        width = width_min if show_var_labels and figwidth < width_min else figwidth
        height = figheight
        
        
    else:
        height_min = len(xy_grouped_data) * 0.2
        height = height_min if show_var_labels and figheight < height_min else figheight
        width = figwidth

    fig, axs = setup_figure(width, height, x_groupby, y_groupby)
    add_xgroup_annotation(x_group_info, axs, x_group_palettes)
    add_ygroup_annotation(y_group_info, axs, y_group_palettes)
    
    if 'legend' in axs:
        patch_dim_ratio = height / width
        set_axes_extents(
            fig,
            axs['legend'],
            x_group_info,
            y_group_info,
            x_group_palettes,
            y_group_palettes,
            # patch_dim_ratio
        )
        add_legend(
            fig,
            axs['legend'],
            x_group_info,
            y_group_info,
            x_group_palettes,
            y_group_palettes,
            # patch_dim_ratio
        )
        
    im = axs['heatmap'].imshow(
        xy_grouped_data,
        cmap = cmap,
        vmin = vmin if vmin else xy_grouped_data.values.min(), 
        vmax = vmax if vmax else xy_grouped_data.values.max(),
        aspect = 'auto',
        interpolation = 'none'
    )
    plt.colorbar(im, axs['colorbar'])
    
    if swap_axes:
        axs['y_axis'].set_yticks([])
    
        if show_var_labels:
            labels = xy_grouped_data.columns
            axs['x_axis'].set_xticks(range(len(labels)))
            axs['x_axis'].set_xticklabels(labels, rotation = 90)

        else:
            axs['x_axis'].set_xticks([])
    
    else:
        axs['x_axis'].set_xticks([])
    
        if show_var_labels:
            labels = xy_grouped_data.index
            axs['y_axis'].set_yticks(range(len(labels)))
            axs['y_axis'].set_yticklabels(labels)

        else:
            axs['y_axis'].set_yticks([])

    return fig, axs, x_group_info, y_group_info
