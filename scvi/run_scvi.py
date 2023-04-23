"""
Train scVI on h5ad containing counts at various hidden dimensions and latent space dimensions.
"""

import scvi
import scanpy as sc
import argparse
import pandas as pd

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
                    prog='run_scvi.py',
                    description='Runs scVI at specified parameters',
                    epilog='Thanks for using!')
    
    parser.add_argument('filename')           # positional argument
    parser.add_argument('-o','--output')
    parser.add_argument('-ld', '--latent', type=int)      # option that takes a value
    parser.add_argument('-hd', '--hidden', type=int)  # on/off flag
    
    args = parser.parse_args()
    
    filename = args.filename
    num_latent = args.latent
    num_hidden = args.hidden
    output = args.output

    adata = sc.read(filename)
    print(adata.shape)

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        continuous_covariate_keys=["pct_counts_mt"],
    )
    
    model = scvi.model.SCVI(
        adata,
        n_hidden=num_hidden,
        n_latent=num_latent
    )
    
    model.train()
    
    latent = model.get_latent_representation()
    model.save(output + "model_{}_{}/".format(num_hidden, num_latent), overwrite=True)
    pd.DataFrame(latent).to_csv(output + 'latent_{}_{}.csv'.format(num_hidden, num_latent))
    
    
    