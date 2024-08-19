import datasets

local_dir='/home/than/Datasets/HEST_bench/' # hest will be dowloaded to this folder

# Note that the full dataset is around 1TB of data

dataset = datasets.load_dataset(
    'MahmoodLab/hest-bench', 
    cache_dir=local_dir,
    # patterns='*'
)