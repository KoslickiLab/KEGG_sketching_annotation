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

2. `dump_hashes.sh` is script that will dump all the hashes of a given sourmash signature collection.
Such as one output from `sketch_genes_KO.sh`. This will be used in the future to create bloom filters or the like
to store the hashes for fast lookups.


## Other sourmash functionality
1. Can use `sourmash sig kmers` along with the training database to extract all the k-mers and the hashes
that were sketched. For example:
```bash
$ sourmash sig kmers --signatures kegg_genes_KO_100.faa_scale_10.db.zip --sequences /data/shared_data/KEGG_data/kegg_genes_KO_100.faa --save-kmers kmer-matches.csv --protein -k 7
```