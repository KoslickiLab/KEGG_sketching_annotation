#!/usr/bin/env python
import argparse
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
# for relative imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from src.HelperFuncs import make_sketches, run_sourmash_gather, check_extension, calculate_sourmash_performance


def main():
    parser = argparse.ArgumentParser(description="This script will use Sourmash gather to classify a metagenome and report statistics on it.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference_file', type=str, help="The reference file (fna or faa) that you want to compare against.")
    parser.add_argument('-m', '--metagenome', type=str, help="The simulated metagenome.")
    parser.add_argument('-o', '--out_dir', type=str, help="The output directory.")
    parser.add_argument('-k', '--kmer_size', type=int, help="The size of the kmer to use.")
    parser.add_argument('-t', '--threshold_bp', type=int, help="The threshold of bp in common between query and "
                                                               "reference to warant inclusion in the results.", default=100)
    #parser.add_argument('-n', '--num_results', type=int, help="Return at most n results", default=1000)
    parser.add_argument('--ref_scale_size', type=int, help="The scale factor to use for the reference database: "
                                                           "s is an integer >=1 and is the denominator of the fraction of sketches to keep.")
    parser.add_argument('--query_scale_size', type=int, help="The scale factor to use for the query: "
                                                           "s is an integer >=1 and is the denominator of the fraction of sketches to keep.")
    parser.add_argument('--query_protein', action='store_true', help="Flag indicating that the query is protein. Otherwise assume DNA.")
    parser.add_argument('--query_translate', action='store_true',
                        help="Flag indicating that the query is DNA, but should be translated to protein.")
    parser.add_argument('--reuse_query_sketch', action='store_true',
                        help="Flag saying the the simulation/query has not changed, so do not update it")
    # parse the args
    args = parser.parse_args()
    query_is_protein = args.query_protein
    query_translate = args.query_translate
    if query_translate and query_is_protein:
        raise Exception("The flags --query_is_protein and --query_translate are mutually exclusive")
    reference_file = args.reference_file
    reuse_query_sketch = args.reuse_query_sketch
    metagenome_file = args.metagenome
    out_dir = args.out_dir
    ksize = args.kmer_size
    ref_scale = args.ref_scale_size
    query_scale = args.query_scale_size
    threshold_bp = args.threshold_bp
    #num_res = args.num_results
    if query_is_protein:
        query_sketch_type = 'protein'
    else:
        query_sketch_type = 'dna'
    if query_translate:
        query_sketch_type = 'translate'
    # check args
    if not exists(out_dir):
        os.makedirs(out_dir)
    if not exists(metagenome_file):
        raise Exception(f"Input metagenome {metagenome_file} does not appear to exist")
    if ref_scale < 1 or query_scale < 1:
        raise Exception(f"Scale sizes must be greater than or equal to one. You provided: {ref_scale} and {query_scale}")
    # check if the reference file sketch with the right params exists
    ref_sketch_file = os.path.join(out_dir, os.path.basename(f"{reference_file}_k_{ksize}_scale_{ref_scale}.sig"))
    if not exists(ref_sketch_file):
        warnings.warn(f"Sketch file {ref_sketch_file} does not exist, creating it now.")
        sketch_type = check_extension(reference_file)
        make_sketches(ksize, ref_scale, reference_file, sketch_type, out_dir, per_record=True)
    # check if the query file sketch with the right params exists
    query_sketch_file = f"{metagenome_file}_k_{ksize}_scale_{query_scale}.sig"
    #if not exists(query_sketch_file):
    #warnings.warn(f"Sketch file {query_sketch_file} does not exist, creating it now.")
    if not reuse_query_sketch:
        make_sketches(ksize, query_scale, metagenome_file, query_sketch_type, out_dir, per_record=False)
    # Then run sourmash gather
    gather_out_file = os.path.join(out_dir, f"{os.path.basename(query_sketch_file)}_{os.path.basename(ref_sketch_file)}_gather.csv")
    if query_translate:
        run_sourmash_gather(query_sketch_file, ref_sketch_file, gather_out_file, 'protein', num_results=None,
                            threshold_bp=threshold_bp, quiet=False)
    else:
        run_sourmash_gather(query_sketch_file, ref_sketch_file, gather_out_file, query_sketch_type, num_results=None,
                            threshold_bp=threshold_bp, quiet=False)
    # And calculate the results
    #stats = calc_binary_stats_sourmash(rel_abund_file, gather_out_file)
    # TODO: use the find_genes_in_sim.py script to calculate the results
    # print(stats)


if __name__ == "__main__":
    main()
