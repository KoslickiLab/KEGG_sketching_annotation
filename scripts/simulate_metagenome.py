#!/usr/bin/env python
import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess


def main():
    parser = argparse.ArgumentParser(description="This script uses BBTools randomreads.sh to simulate reads from a"
                                                 "metagenome.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference_file', type=str, help="Database that you want to generate reads from.")
    parser.add_argument('-o', '--out_file', type=str, help="The output simulated metagenome.")
    parser.add_argument('-n', '--num_reads', type=int, help="The number of reads to simulate.")
    parser.add_argument('-l', '--len_reads', type=int, help="The length of the reads to generate.")
    parser.add_argument("--noisy", action='store_true', help="The full path to the directory that the files will be written")
    args = parser.parse_args()
    simple_names = "t"
    overwrite = "t"
    illumina_names = "t"


if __name__ == "__main__":
    main()

