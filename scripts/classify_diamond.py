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
from src.HelperFuncs import make_sketches, run_sourmash_gather, check_extension, \
     build_diamond_db, run_diamond_blastx


def main():
    parser = argparse.ArgumentParser(description="This script will use Sourmash gather to classify a metagenome and report statistics on it.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference_file', type=str, help="The reference file (fna or faa) that you want to compare against.")
    parser.add_argument('-m', '--metagenome', type=str, help="The simulated metagenome.")
    parser.add_argument('-o', '--out_dir', type=str, help="The output directory.")
    parser.add_argument('-d', '--diamond_file', type=str, help="The output diamond file.")
    # parse the args
    args = parser.parse_args()
    reference_file = args.reference_file
    metagenome_file = args.metagenome
    out_dir = args.out_dir
    diamond_file = args.diamond_file
    # check args
    if not exists(out_dir):
        os.makedirs(out_dir)
    if not exists(metagenome_file):
        raise Exception(f"Input metagenome {metagenome_file} does not appear to exist")
    # Check if the reference database has been built, and build it if it hasn't
    ref_db = os.path.join(out_dir, os.path.basename(f"{reference_file}.dmnd"))
    if not exists(ref_db):
        warnings.warn(f"Diamond database {ref_db} does not appear to exist. Making it now. This may take some time.")
        build_diamond_db(reference_file, ref_db)
    # Do the alignment
    out_file = diamond_file
    if out_file is None:
        out_file = os.path.join(out_dir, f"{os.path.basename(metagenome_file)}_{os.path.basename(ref_db)}_matches.csv")
    run_diamond_blastx(metagenome_file, ref_db, out_file)
    # Calculate the results
    # TODO: use the find_genes_in_sim.py script to calculate the results
    #stats = calc_binary_stats_diamond(metagenome_file, out_file)
    # Print the results
    #print(stats)

if __name__ == "__main__":
    # Example of running it as a script from the script dir:
    # ./classify_and_report_diamond.py -r ../test_data/input/kegg_genes_KO.faa -m ../test_data/output/test_simulation.fq -o ../test_data/output/
    main()
