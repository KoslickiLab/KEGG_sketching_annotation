# KEGG_sketching_annotation
Scripts to sketch KEGG and explore using FracMinHash as a way to functionally annotate a metagenome.

## Requirements
The `requirements.txt` file contains the necessary packages required to run the code in this repo.
You can install it via:
```commandline
conda create -y --name KEGG_env python=3.8
conda install -y --name KEGG_env -c conda-forge -c bioconda --file requirements.txt
conda activate KEGG_env
```

## Data
The data this repo uses is currently on the lab GPU server, under `/data/shared_data/`.
This was extracted with [this repo](https://github.com/KoslickiLab/KEGG_data_extraction).
