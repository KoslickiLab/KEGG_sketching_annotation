# KEGG scripts
This directory will contain scripts specific to the KEGG data itself

## Data location
The data is located in the `/data/shared_data/KEGG_data` directory on the lab GPU server (e5-cse-cbdmk02).

## Scripts
1. `sketch_genes_KO.sh` is a bash script that will create sourmash sketches for the KEGG genes KO.fna/faa files.
A range of k-mer sizes and scaled values are computed. Namely:

| k-size       | scaled size      | molecule type |
|--------------|------------------|---------------|
| 21, 31, 51   | 1, 10, 100, 1000 | NT            |
| 5, 7, 11, 15 | 1, 10, 100, 1000 | AA            |

2. `dump_hashes.py` is a python script that will dump all the hashes of a given sourmash signature collection.
Such as one output from `sketch_genes_KO.sh`. This will be used in the future to create bloom filters or the like
to store the hashes for fast lookups.