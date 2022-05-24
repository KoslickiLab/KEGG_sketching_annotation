#!/usr/bin/env python
import argparse
import os
from os import listdir
from os.path import isfile, join


def main():
    parser = argparse.ArgumentParser(description="This script creates training/reference sketches of KEGG. "
                                                 "Sketches will be made for each record contained in each *.fna and"
                                                 "*.faa file contained in the --in_dir",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--scale_factor', type=int, help="Denominator of the scale factor to use", default=500)
    parser.add_argument('-k', '--k_size', type=int, help="K-mer size", default=21)
    parser.add_argument("--in_dir", type=str, help="The full path to the directory that contains the FASTA files to sketch")
    parser.add_argument("--out_dir", type=str, help="The full path to the directory that the files will be written")
    args = parser.parse_args()
    ksize = args.k_size
    scale_factor = args.scale_factor
    in_dir = args.in_dir
    out_dir = args.out_dir
    if scale_factor < 1:
        raise Exception(f"Scale factor must be an integer greater than or equal to 1")
    if not os.path.exists(in_dir):
        raise Exception(f"The supplied input directory {in_dir} does not exist.")
    file_names = [os.path.join(in_dir, f) for f in listdir(in_dir) if isfile(join(in_dir, f))]
    if not file_names:
        raise Exception(f"The directory {in_dir} is empty")


if __name__ == "__main__":
    main()

