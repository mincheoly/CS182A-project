# CS182A-project
Using deep latent representations for improving linear hypothesis testing

# Steps to reproduce the results

1. Follow the notebook preprocess.ipynb to generate AnnData objects from raw data.
2. Use submit_scvi.sh in the scvi folder to generate the latent variables. 
3. Follow the setup_analysis.ipynb notebook in the analysis folder to setup the paired guide analysis.
4. Use run_regressions.sh in the analysis folder to compute effect sizes.
5. Follow the paired_guide.ipynb notebook in the analysis folder to plot and investigate the positive control metrics.
