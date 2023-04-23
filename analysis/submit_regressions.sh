# python run_regressions.py \
#     --h5ad_path /mnt/c/Users/minch/Data/bmdc/processed_bmdc_counts.h5ad \
#     --design_path /mnt/c/Users/minch/Data/bmdc/analysis/design_matrix.csv \
#     --output /mnt/c/Users/minch/Data/bmdc/analysis/beta_null.csv

python run_regressions.py \
    --h5ad_path /mnt/c/Users/minch/Data/bmdc/processed_bmdc_counts.h5ad \
    --latent_path /mnt/c/Users/minch/Data/bmdc/analysis/latent_128_5.csv \
    --design_path /mnt/c/Users/minch/Data/bmdc/analysis/design_matrix.csv \
    --output /mnt/c/Users/minch/Data/bmdc/analysis/beta_128_5.csv

# python run_regressions.py \
#     --h5ad_path /mnt/c/Users/minch/Data/bmdc/processed_bmdc_counts.h5ad \
#     --latent_path /mnt/c/Users/minch/Data/bmdc/analysis/latent_0_0.csv \
#     --design_path /mnt/c/Users/minch/Data/bmdc/analysis/design_matrix.csv \
#     --output /mnt/c/Users/minch/Data/bmdc/analysis/beta_0_0.csv