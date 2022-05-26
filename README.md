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

## Main idea
Sketch the AA and NT sequences, compare the following approaches:

1. AA (vary the k-mer size and sketch size) with `sourmash gather`
2. NT (vary the k-mer size and sketch size) with `sourmash gather`
3. AA/NT prefetch, then alignment
4. Cmash/mash screen approach: sketch the database, stream the query over it looking for all matches

Using the following query data (in increasing level of realism):
1. Simulated subsequences of the AA/NT database
2. Simulated metagenomes that have AA/NT database sequences spiked in
3. Full simulated metagenomes (AA+prodigal/AA+6-frame-translation/NT)
4. Mock metagenomes
5. Real metagenomes that have been extensively analyzed previously

## Structure of repo
`src` contains the module with various helper functions, wrappers to sourmash, etc.
`scripts` contains the scripts that perform various tasks (sketching, simulation, etc.)
`utils` contains third party tools (in this case, bbtools) for things like making simulations.
`test_data` contains a small portion of real KEGG data for rapid testing and development.

## Order of scripts to run

1. `make_sketches.py` This will take the FAA and FNA files and create sketches/signatures
for each one.
2. `simulate_metagenome.py` simulates a metagenome when provided with a reference database. This will need some work since currently bbtools uses the _whole_ database to simulate a metagenome (i.e. coverage is very low), so we'll need to do some down-sampling. This will also create a ground-truth functional profile (though that needs to be checked for accuracy).
3. `classify_and_report.py` is intended to run sourmash gather, parse the results, and then compare this to the ground truth.
