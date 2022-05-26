#!/usr/bin/env python
import argparse
import os
import sys
from os import listdir
from os.path import isfile, join
# for relative imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from src.HelperFuncs import make_sketches, check_extension


def main():
    parser = argparse.ArgumentParser(description="This script creates training/reference sketches of KEGG. "
                                                 "Sketches will be made for each record contained in each *.fna and"
                                                 "*.faa file contained in the --in_dir. If the --in_dir is not specified,"
                                                 "the single file --single_file will be sketched instead.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--scale_factor', type=int, help="Denominator of the scale factor to use", default=500)
    parser.add_argument('-k', '--k_size', type=int, help="K-mer size", default=21)
    parser.add_argument("--in_dir", type=str, help="The full path to the directory that contains the FASTA files to sketch")
    parser.add_argument("--out_dir", type=str, help="The full path to the directory that the files will be written")
    parser.add_argument("--single_file", type=str, help="The full path to a single file you wish to sketch")
    parser.add_argument("--protein", action='store_true', help="Flag if the --single_file you provide is a protein file.")
    args = parser.parse_args()
    ksize = args.k_size
    scale_factor = args.scale_factor
    in_dir = args.in_dir
    out_dir = args.out_dir
    single_file = args.single_file
    if single_file:
        if not os.path.exists(single_file):
            raise Exception(f"The file {single_file} does not exist")
        if single_file and in_dir:
            raise Exception("The parameters --in_dir and --single_file are mutually exclusive. Please specify just one")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if scale_factor < 1:
        raise Exception(f"Scale factor must be an integer greater than or equal to 1")
    # If you specified a directory, sketch everything in it
    if in_dir:
        if not os.path.exists(in_dir):
            raise Exception(f"The supplied input directory {in_dir} does not exist.")
        file_names = [os.path.abspath(os.path.join(in_dir, f)) for f in listdir(in_dir) if isfile(join(in_dir, f))]
        if not file_names:
            raise Exception(f"The directory {in_dir} is empty")
        for file_name in file_names:
            prefix, suffix = file_name.strip().split(".")
            if suffix not in ["fna", "faa"]:
                raise Exception(f"The file {file_name} does not have the suffix faa, fna. Cannot proceed as I don't "
                                f"know if this is a nucleotide sequence or amino acid sequence.")
        # Now actually do the sketching
        for file_name in file_names:
            sketch_type = check_extension(file_name)
            print(f"Sketching the entries in the file {file_name}")
            make_sketches(ksize, scale_factor, file_name, sketch_type, out_dir, per_record=True)
    # Otherwise you specified a single file, so sketch just that one
    elif single_file:
        file_name = single_file
        if args.protein:
            sketch_type = "aa"
        else:
            sketch_type = "nt"
        print(f"Sketching the entries in the file {file_name}")
        make_sketches(ksize, scale_factor, file_name, sketch_type, out_dir, per_record=False)


if __name__ == "__main__":
    main()

