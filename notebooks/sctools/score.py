import numpy as np
import pandas as pd


def bin_data(average_expression, nbins):
    n_values = len(average_expression)
    num_vals_per_bin = n_values // nbins
    sorted_idx = np.argsort(average_expression)
    
    bin_assignments = np.zeros(shape = n_values, dtype = int)
    for bin_id, cut_idx in enumerate(range(0, n_values, num_vals_per_bin)):
        value_idx = sorted_idx[cut_idx: cut_idx + num_vals_per_bin]
        
        if len(value_idx) < num_vals_per_bin:
            bin_assignments[value_idx] = bin_id - 1
            break
            
        bin_assignments[value_idx] = bin_id
    
    return bin_assignments


# gene module score as defined by Tirosh et al. Science 2016
# https://doi.org/10.1126/science.aad0501
def gene_module_score(adata, gene_list):
    genes_in_adata = set(adata.var.index.to_list())
    filtered_gene_list = list(
        set(gene_list) & genes_in_adata
    )
    
    average_expression = np.array(adata.X.mean(axis = 0)).flatten()
    gene_bins = pd.Series(
        bin_data(
            average_expression, 
            25
        ),
        index = adata.var.index
    )
    
    control_genes = set()
    for gene in filtered_gene_list:
        bin_id = gene_bins[gene]
        genes_in_bin = gene_bins[gene_bins == bin_id].index
        random_control_genes = np.random.choice(
            genes_in_bin, 
            size = 100, 
            replace = False
        )
        control_genes.update(random_control_genes)
    
    control_scores = np.array(adata[:, list(control_genes)].X.mean(axis = 1)).flatten()
    gene_list_scores = np.array(adata[:, filtered_gene_list].X.mean(axis = 1)).flatten()
    
    return gene_list_scores - control_scores 