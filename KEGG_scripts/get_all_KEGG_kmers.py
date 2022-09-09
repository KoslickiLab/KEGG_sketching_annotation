#!/usr/bin/env python
import argparse
import screed


# parse arguments
parser = argparse.ArgumentParser(description='Get all kmers from a fasta file. This is intended for use with KEGG amino '
                                             'acid sequences and for relatively small k-mer sizes (e.g. 3-6).')
parser.add_argument('-f', '--fasta', help='input KEGG FAA fasta file')
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
        seq_id = seq_name.split('|')[0]
        seq_kegg = seq_name.split('|')[2]  #TODO: not using this atm, since can easily map id <-> KO
        # get the k-mers
        kmers = set(seq[i:i + ksize] for i in range(len(seq) - ksize + 1))
        # add the k-mers to the dictionary
        for kmer in kmers:
            if kmer not in kmer_dict:
                kmer_dict[kmer] = set([seq_id])
            else:
                kmer_dict[kmer].add(seq_id)

# write the k-mers to a file
with open(out_file, 'w') as f:
    for kmer in kmer_dict:
        f.write(f"{kmer}, {'|'.join(kmer_dict[kmer])}\n")
