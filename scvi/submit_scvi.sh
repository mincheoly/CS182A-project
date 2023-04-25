python run_scvi.py \
    /mnt/c/Users/minch/Data/bmdc/processed_bmdc_hvg_counts.h5ad \
    -o /mnt/c/Users/minch/Data/bmdc/analysis/ \
    -ld 5 \
    -hd 128
    
python3 scvi/run_scvi.py ./processed_bmdc_hvg_counts.h5ad -o analysis/ -ld 5 -hd 64
python3 scvi/run_scvi.py ./processed_bmdc_hvg_counts.h5ad -o analysis/ -ld 10 -hd 64
python3 scvi/run_scvi.py ./processed_bmdc_hvg_counts.h5ad -o analysis/ -ld 20 -hd 64
python3 scvi/run_scvi.py ./processed_bmdc_hvg_counts.h5ad -o analysis/ -ld 5 -hd 128
python3 scvi/run_scvi.py ./processed_bmdc_hvg_counts.h5ad -o analysis/ -ld 10 -hd 128
python3 scvi/run_scvi.py ./processed_bmdc_hvg_counts.h5ad -o analysis/ -ld 20 -hd 128