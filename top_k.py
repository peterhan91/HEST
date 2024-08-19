import os 
import json
import pandas as pd
import numpy as np
import scanpy as sc
from typing import List
from tqdm import tqdm

def get_k_genes_from_df(meta_df: pd.DataFrame, k: int, criteria: str, save_dir: str=None) -> List[str]:
    """Get the k genes according to some criteria across common genes in all the samples in the HEST meta dataframe

    Args:
        meta_df (pd.DataFrame): HEST meta dataframe
        k (int): number of genes to return
        criteria (str): criteria for the k genes to return
            - 'mean': return the k genes with the largest mean expression across samples
            - 'var': return the k genes with the largest expression variance across samples
        save_dir (str, optional): genes are saved as json array to this path if not None. Defaults to None.

    Returns:
        List[str]: k genes according to the criteria
    """
    
    adata_list = []
    for _, row in meta_df.iterrows():
        id = row['id']
        adata = sc.read_h5ad(f"/home/than/Datasets/HEST_data/st/{id}.h5ad")
        adata_list.append(adata)
    return get_k_genes(adata_list, k, criteria, save_dir=save_dir)


def get_k_genes(adata_list: List[sc.AnnData], k: int, criteria: str, save_dir: str=None, min_cells_pct=0.10) -> List[str]: # type: ignore
    """Get the k genes according to some criteria across common genes in all the samples in the adata list

    Args:
        adata_list (List[sc.AnnData]): list of scanpy AnnData containing gene expressions in adata.to_df()
        k (int): number of most genes to return
        criteria (str): criteria for the k genes to return
            - 'mean': return the k genes with the largest mean expression across samples
            - 'var': return the k genes with the largest expression variance across samples
        save_dir (str, optional): genes are saved as json array to this path if not None. Defaults to None.
        min_cells_pct (float): filter out genes that are expressed in less than min_cells_pct% of the spots for each slide

    Returns:
        List[str]: k genes according to the criteria
    """
    # check_arg(criteria, 'criteria', ['mean', 'var'])
    
    common_genes = None
    stacked_expressions = None

    # Get the common genes
    for adata in adata_list:
        my_adata = adata.copy()
        
        if min_cells_pct:
            print('min_cells is ', np.ceil(min_cells_pct * len(my_adata.obs)))
            sc.pp.filter_genes(my_adata, min_cells=np.ceil(min_cells_pct * len(my_adata.obs)))
        curr_genes = np.array(my_adata.to_df().columns)
        if common_genes is None:
            common_genes = curr_genes
        else:
            common_genes = np.intersect1d(common_genes, curr_genes)
            

    common_genes = [gene for gene in common_genes if 'BLANK' not in gene and 'Control' not in gene]

    for adata in tqdm(adata_list):

        if stacked_expressions is None:
            stacked_expressions = adata.to_df()[common_genes]
        else:
            stacked_expressions = pd.concat([stacked_expressions, adata.to_df()[common_genes]])
    
    stacked_expressions.to_csv('stacked_expressions.csv')
    
    if criteria == 'mean':
        nb_spots = len(stacked_expressions)
        mean_expression = stacked_expressions.sum() / nb_spots
        
        top_k = mean_expression.nlargest(k).index
    elif criteria == 'var':
        # stacked_adata = sc.AnnData(stacked_expressions.astype(np.float32))
        stacked_adata = sc.AnnData(stacked_expressions)
        sc.pp.filter_genes(stacked_adata, min_cells=0)
        sc.pp.log1p(stacked_adata)
        sc.pp.highly_variable_genes(stacked_adata, n_top_genes=k)
        top_k = stacked_adata.var_names[stacked_adata.var['highly_variable']][:k].tolist()
    else:
        raise NotImplementedError()

    if save_dir is not None:
        json_dict = {'genes': list(top_k)}
        with open(save_dir, mode='w') as json_file:
            json.dump(json_dict, json_file)

    print(f'selected genes {top_k}')
    return top_k



if __name__ == '__main__':
    meta_df = pd.read_csv("hf://datasets/MahmoodLab/hest/HEST_v1_0_2.csv")
    get_k_genes_from_df(meta_df, 200, 'var', 'top_genes.json')