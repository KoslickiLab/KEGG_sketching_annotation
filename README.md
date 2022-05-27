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


## Classification results
Running the following will print out binary classification results:
```commandline
./classify_and_report.py -r ../test_data/input/kegg_genes_KO.fna -m ../test_data/output/test_simulation.fq -o ../test_data/output/ -k 21 -t 100 -n 1000 --ref_scale_size 100 --query_scale_size 100
```
which will result in something like:
```commandline
{'TP': 97, 'FP': 0, 'FN': 3, 'precision': 1.0, 'recall': 0.97, 'F1': 0.9847715736040609}
```

# Comparison to DIAMOND
DIAMOND is an alignment-based functional annotation software package that is widely used.
Let's compare the sourmash/FracMinHash approach to DIAMOND. 

## Install instructions
Updated from [here](https://github.com/bbuchfink/diamond/wiki):
```commandline
# downloading the tool
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
# creating a diamond-formatted database file
./diamond makedb --in reference.fasta -d reference
# running a search in blastp mode
./diamond blastp -d reference -q queries.fasta -o matches.tsv
# running a search in blastx mode
./diamond blastx -d reference -q reads.fasta -o matches.tsv
# downloading and using a BLAST database
update_blastdb.pl --decompress --blastdb_version 5 swissprot
./diamond prepdb -d swissprot
./diamond blastp -d swissprot -q queries.fasta -o matches.tsv
```

## Running DIAMOND and parsing the results
After making the simulation, run via:
```commandline
 ./classify_and_report_diamond.py -r ../test_data/input/kegg_genes_KO.faa -m ../test_data/output/test_simulation.fq -o ../test_data/output/
```
And you should get a result like
```commandline
{'TP': 98, 'FP': 110, 'FN': 2, 'precision': 0.47115384615384615, 'recall': 0.98, 'F1': 0.6363636363636364, 'Percent correct alignments': 0.8630737190242755, 'Total number of alignments': 1156228, 'Total number of sequences': 1000000.0}
```

# CAMISIM
Following the guide [here](https://github.com/KoslickiLab/useful_tools/tree/main/Metagenomics/simulate_metagenomic_by_CAMISIM) (thanks Shaopeng!) to set things up.

# Other methods to compare against
https://github.com/biobakery/shortbred
