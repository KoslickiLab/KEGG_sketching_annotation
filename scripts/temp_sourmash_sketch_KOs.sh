#!/usr/bin/env bash
set -e
set -u
# where your input data is
DATADIR=/data/shared_data/Classifying_to_KEGG/data/metagenomes/
# where you want the output written
OUTDIR=/data/shared_data/Classifying_to_KEGG/output/
# the metagenome to bo classified
FILEPREFIX="SRS1041031.denovo_duplicates_marked.trimmed"
# the desired file name for the ouput of the sketched metagenome
SIGFILE=${DATADIR}/${FILEPREFIX}.sig.zip
# the reference database containing the sketched KOs
REF=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_10.sig.zip
# sketch the metagenome
sourmash sketch translate -f -p scaled=1000,k=5 ${DATADIR}/${FILEPREFIX}*.fastq -o ${SIGFILE} --merge ${SIGFILE}
# infer KO frequencies. Results will be in the *gather_k_5.csv file in your specified OUTDIR
sourmash gather --protein -k 5 --estimate-ani-ci --threshold-bp 500 ${SIGFILE} ${REF} -o ${OUTDIR}/${FILEPREFIX}_$(basename $REF)_gather_k_5.csv
