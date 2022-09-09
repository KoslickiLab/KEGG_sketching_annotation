#!/usr/bin/env python
import argparse
import screed
import os
import sys
from os.path import exists
import pathlib
from os import listdir
from os.path import isfile, join
import subprocess
import re
import numpy as np
from collections import Counter
import warnings
import pandas as pd

# parse arguments
parser = argparse.ArgumentParser(description='Get all kmers from a fasta file. This is intended for use with amino '
                                             'acid sequences and for relatively small k-mer sizes (e.g. 3-6).')
parser.add_argument('-f', '--fasta', help='input fasta file')
parser.add_argument('-k', '--ksize', help='kmer size', default=5)
parser.add_argument('-o', '--out', help='output file. Will be a plain text file, csv separated by newlines.')
args = parser.parse_args()
ksize = int(args.ksize)
fasta_file = args.fasta
out_file = args.out
# make sure ksize is small
if ksize > 7:
    raise ValueError('ksize must be <= 7')

# initialize the kmer counter, keys are k-mers, values are sequence names
kmer_dict = {}
# open the fasta file
with screed.open(fasta_file) as seqfile:
    for record in seqfile:
        # get the sequence
        seq = record.sequence
        # get the sequence name
        seq_name = record.name
        # For the sequences from the KEGG sketching repo, the first token is the sequence name
        seq_name = seq_name.split('|')[0]
        # get the k-mers
        kmers = set(seq[i:i + ksize] for i in range(len(seq) - ksize + 1))
        # add the k-mers to the dictionary
        for kmer in kmers:
            if kmer not in kmer_dict:
                kmer_dict[kmer] = set(seq_name)
            else:
                kmer_dict[kmer].add(seq_name)

# write the k-mers to a file
with open(out_file, 'w') as f:
    for kmer in kmer_dict:
        print(f"header: {'|'.join(kmer_dict[kmer])}")
        print(f"val: {kmer_dict[kmer]}")
        f.write(f"{kmer}, {'|'.join(kmer_dict[kmer])}\n")
