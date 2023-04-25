"""
Get reconstruction error of saved scVI models using h5ad containing counts.
"""

import scvi
import scanpy as sc
import argparse
import pandas as pd
import numpy as np

if _name_ == '_main_':
    
    parser = argparse.ArgumentParser(
                    prog='run_scvi.py',
                    description='Gets reconstruction error of scVI models at specified parameters',
                    epilog='Thanks for using!')
    
    parser.add_argument('filename')                       # h5ad file filename
    parser.add_argument('-o','--output')
    parser.add_argument('-ld', '--latent', nargs="+", type=int)      # option that takes a list of values for latent dimension
    parser.add_argument('-hd', '--hidden', nargs="+", type=int)      # option that takes a list of values for hidden dimension
    
    args = parser.parse_args()
    
    filename = args.filename
    latent = args.latent
    hidden = args.hidden
    output = args.output

    adata = sc.read(filename)

    reconstruction_errors = np.zeros((len(hidden), len(latent)))

    for i, num_hidden in enumerate(hidden):
        for j, num_latent in enumerate(latent):
            model = scvi.model.SCVI.load("analysis/model_{}_{}/".format(num_hidden, num_latent), adata = adata)
            r = model.get_reconstruction_error()
            reconstruction_errors[i][j] = model.get_reconstruction_error().get("reconstruction_loss")
    df = pd.DataFrame(reconstruction_errors, index = hidden, columns = latent)
    df.to_csv(output + 'reconstruction_errors.csv')