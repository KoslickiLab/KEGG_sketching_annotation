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
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from src.HelperFuncs import make_sketches, run_sourmash_gather, check_extension, calculate_sourmash_performance


def main():
    parser = argparse.ArgumentParser(description="This script will use return performance statistics when comparing sourmash to the ground truth.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sourmash_results', type=str, help="The output csv from sourmash gather")
    parser.add_argument('-g', '--ground_truth', type=str, help="The ground truth as calculated by find_genes_in_sim.py.")
    parser.add_argument('-o', '--out_file', type=str, help="The output csv file containing the performance statistics.")
    # parse the args
    args = parser.parse_args()
    sourmash_results = args.sourmash_results
    ground_truth = args.ground_truth
    out_file = args.out_file
    # check if files exist
    if not exists(sourmash_results):
        raise Exception(f"The sourmash results file does not exist: {sourmash_results}")
    if not exists(ground_truth):
        raise Exception(f"The ground truth file does not exist: {ground_truth}")
    # get the performance statistics
    performance_df = pd.DataFrame()
    for filter_threshold in range(0, 100, 5):
        # append the row to the dataframe
        performance_df = pd.concat([performance_df, calculate_sourmash_performance(sourmash_results, ground_truth, filter_threshold)], ignore_index=True)

    # write the performance statistics to a csv file
    performance_df.to_csv(out_file, index=False)



if __name__ == "__main__":
    main()
