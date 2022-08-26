#!/usr/bin/env bash
set -e
set -u
KEGGDATAAA=/data/shared_data/KEGG_data/kegg_genes_KO.faa
KEGGDATANT=/data/shared_data/KEGG_data/kegg_genes_KO.fna
OUTDIR=/data/shared_data/KEGG_data/sourmash_sketches/output/
for sketchSize in 1 10 100 1000
do
        /usr/bin/time sourmash sketch dna -f -o ${OUTDIR}/$(basename $KEGGDATANT)_scale_${sketchSize}.db.zip --singleton -p k=21,k=31,k=51,abund,scaled=${sketchSize} $KEGGDATANT &
        /usr/bin/time sourmash sketch protein -f -o ${OUTDIR}/$(basename $KEGGDATAAA)_scale_${sketchSize}.db.zip --singleton -p k=5,k=7,k=11,k=15,abund,scaled=${sketchSize} $KEGGDATAAA &
done
