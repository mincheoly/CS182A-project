import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import scipy.sparse as sparse
import scipy.stats as stats
import sklearn.linear_model as lm
import sklearn.metrics as metrics
import itertools
import seaborn as sns


def normalize(adata):
    
    adata.layers['counts'] = adata.X
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    

def corr2_coeff(A, B):
    """
    From: https://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays
    """
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) /ssA[:, None]


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
                prog='run_scvi.py',
                description='Runs scVI at specified parameters',
                epilog='Thanks for using!')
    
    parser.add_argument('-h5ad', '--h5ad_path')           # positional argument
    parser.add_argument('-latent', '--latent_path')
    parser.add_argument('-design', '--design_path')
    parser.add_argument('-o','--output')

    args = parser.parse_args()
    
    h5ad_path = args.h5ad_path
    latent_path = args.latent_path
    design_path = args.design_path
    output = args.output
    
    # Read and normalize anndata
    adata = sc.read(h5ad_path)
    adata = adata[adata.obs['tp'] == '3hr'].copy()

    normalize(adata)
    num_genes = adata.shape[1]
    
    # Read the cell state matrix
    cell_state_matrix = pd.read_csv(latent_path).values
    design_matrix = pd.read_csv(design_path, index_col=0)
    guide_list = design_matrix.columns.tolist()
    num_guides = len(guide_list)
    
    # Converting all genes to dense is too memory intensive
    num_genes_in_chunk = 1000
    num_states = 5
    num_chunks = int(num_genes/num_genes_in_chunk)+1
    expr = adata.X
    
    coefs = np.zeros((num_genes, num_guides))
    for chunk in range(num_chunks):
        print('Working on chunk', chunk, 'of', num_chunks)
        start = num_genes_in_chunk*chunk
        end = start+num_genes_in_chunk
        y = expr[:, start:end].toarray()
        cell_state_model = lm.LinearRegression(n_jobs=5).fit(cell_state_matrix,y)
        y_partial = y-cell_state_model.predict(cell_state_matrix)

        beta = corr2_coeff(design_matrix.values.T, y_partial.T).T

        coefs[start:end, :] = beta
    
    coefs_df = pd.DataFrame(coefs, columns=guide_list, index=adata.var.index)
    
    coefs_df.to_csv(output)