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
    # parse the args
    args = parser.parse_args()
    reference_file = args.reference_file
    metagenome_file = args.metagenome
    out_dir = args.out_dir
    ksize = args.kmer_size
    ref_scale = args.ref_scale_size
    query_scale = args.query_scale_size
    # check args
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(metagenome_file):
        raise Exception(f"Input metagenome {metagenome_file} does not appear to exist")
    if ref_scale < 1 or query_scale < 1:
        raise Exception(f"Scale sizes must be greater than or equal to one. You provided: {ref_scale} and {query_scale}")

    # TODO: check if the reference file has been sketched to the appropriate size, if not, do it on the fly. Same for the metagenome
    # TODO: check if the relative abundances have already been computed, if not, do them on the fly
    # TODO: then run sourmash gather, save the results, and compute binary metrics using the results





if __name__ == "__main__":
    main()
