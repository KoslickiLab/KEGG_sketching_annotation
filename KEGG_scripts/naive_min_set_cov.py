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

# This is a simple min-set-cov implementation for use with short amino acid sequences.

def get_all_kmers(fasta_file, ksize):
    """
    Get all kmers from a fasta file. This is intended for use with amino acid sequences and for relatively small k-mer sizes (e.g. 3-6).
    :param fasta_file: input fasta file
    :param ksize: kmer size
    """
    # make sure ksize is small
    if ksize > 7:
        raise ValueError('ksize must be <= 7')

    # initialize the kmer counter, keys are k-mers, values are sequence names
    kmer_set = set()
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
            kmer_set.update(kmers)

    return kmer_set

intersection_threshold = 10


# parse arguments
parser = argparse.ArgumentParser(description='Run a simple min-set-cov implementation for use with short amino acid sequences.')
parser.add_argument('-f', '--fasta', help='input fasta query/metagenome file')
parser.add_argument('-r', '--ref', help='reference database')
parser.add_argument('-o', '--out', help='output file. Will be a plain text file, csv separated by newlines.')

args = parser.parse_args()
fasta_file = args.fasta
ref_file = args.ref
out_file = args.out

# import the reference database
kmer_to_seq = dict()
with open(ref_file, 'r') as f:
    for line in f.readlines():
        kmer, seqs = line.strip().split(',')
        kmer_to_seq[kmer] = set(seqs.split('|'))

# Reverse this dictionary so that the keys are the sequences and the values are the k-mers
seq_to_kmer = dict()
for kmer in kmer_to_seq:
    for seq in kmer_to_seq[kmer]:
        if seq not in seq_to_kmer:
            seq_to_kmer[seq] = {kmer}
        else:
            seq_to_kmer[seq].add(kmer)

# get the all the kmers from the query
query_kmers = get_all_kmers(fasta_file, 5)

# loop through all the k-mers in the query and update the sequences that have been seen
seq_2_count = dict()
for kmer in query_kmers:
    if kmer in kmer_to_seq:
        for seq in kmer_to_seq[kmer]:
            if seq not in seq_2_count:
                seq_2_count[seq] = 1
            else:
                seq_2_count[seq] += 1

# take the largest seq found
max_seq = max(seq_2_count, key=seq_2_count.get)
max_count = seq_2_count[max_seq]
while max_count >= intersection_threshold:
    print(f"Max seq: {max_seq} with {seq_2_count[max_seq]} kmers in common with the query.")
    # Get the k-mers that are in common with the query and this max seq
    kmer_overlap = seq_to_kmer[max_seq].intersection(query_kmers)
    # remove the k-mers that are in common with the query and this max seq
    query_kmers = query_kmers.difference(kmer_overlap)
    # set the current max seq count to zero
    seq_2_count[max_seq] = 0
    # then repeat
    max_seq = max(seq_2_count, key=seq_2_count.get)
    max_count = seq_2_count[max_seq]
