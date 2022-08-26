#!/usr/bin/env bash
set -e
set -u

# This is for the small test data
KEGGDATAAA=/data/shared_data/KEGG_data/kegg_genes_KO_100.faa
KEGGDATANT=/data/shared_data/KEGG_data/kegg_genes_KO_100.fna
OUTDIR=/data/shared_data/KEGG_data/sourmash_sketches/output/test
SKETCH=/data/shared_data/KEGG_data/sourmash_sketches/output/test/kegg_genes_KO_100.faa_scale_10.db.zip

sourmash sig kmers --signatures $SKETCH --sequences KEGGDATAAA --save-kmers kmer-matches.csv --protein -k 7