#!/usr/bin/env python
import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess
import re
import numpy as np
from collections import Counter
from ..src.HelperFuncs import make_sketches, run_simulation, compute_rel_abundance

def main():
    parser = argparse.ArgumentParser(description="This script will use Sourmash gather to classify a metagenome and report statistics on it.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference_file', type=str, help="The reference file (fna or faa) that you want to compare against.")
    parser.add_argument('-m', '--metagenome', type=str, help="The simulated metagenome.")
    parser.add_argument('-o', '--out_dir', type=str, help="The output directory.")
    parser.add_argument('-k', '--kmer_size', type=int, help="The size of the kmer to use.")
    parser.add_argument('--ref_scale_size', type=int, help="The scale factor to use for the reference database: "
                                                           "s is an integer >=1 and is the denominator of the fraction of sketches to keep.")
    parser.add_argument('--query_scale_size', type=int, help="The scale factor to use for the query: "
                                                           "s is an integer >=1 and is the denominator of the fraction of sketches to keep.")
    parser.add_argument("--noisy", action='store_true', help="The full path to the directory that the files will be written")
    # parse the args
    args = parser.parse_args()
    reference_file = args.reference_file
    out_dir = args.out_dir



if __name__ == "__main__":
    main()
