#!/usr/bin/env python
import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess


def make_sketches(ksize, scale_factor, file_name, sketch_type, out_dir, per_record=False):
    """
    This helper function will create the signature/sketches using sourmash
    :param ksize: the k-size to use
    :param scale_factor: the denominator of the scale factor to use (so >=1)
    :param file_name: the file to sketch
    :param sketch_type: amino acid (aa) or nucleotide (nt)
    :param out_dir: Where to write the signature
    :param per_record: If you want sketches of each entry in the fasta file, or just of the full fasta file (default: False)
    :return: None
    """
    out_file = f"{os.path.join(out_dir, os.path.basename(file_name))}_k_{ksize}_scale_{scale_factor}.sig"
    if per_record:
        cmd = f"sourmash sketch {sketch_type} -p k={ksize},scaled={scale_factor},abund -o {out_file} --singleton {file_name}"
    else:
        cmd = f"sourmash sketch {sketch_type} -p k={ksize},scaled={scale_factor},abund -o {out_file} --singleton {file_name}"
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)


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
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if scale_factor < 1:
        raise Exception(f"Scale factor must be an integer greater than or equal to 1")
    if not os.path.exists(in_dir):
        raise Exception(f"The supplied input directory {in_dir} does not exist.")
    file_names = [os.path.abspath(os.path.join(in_dir, f)) for f in listdir(in_dir) if isfile(join(in_dir, f))]
    if not file_names:
        raise Exception(f"The directory {in_dir} is empty")
    for file_name in file_names:
        prefix, suffix = file_name.strip().split(".")
        if suffix not in ["fna", "faa"]:
            raise Exception(f"The file {file_name} does not have the suffix faa or fna. Cannot proceed as I don't "
                            f"know if this is a nucleotide sequence or amino acid sequence.")
    # Now actually do the sketching
    for file_name in file_names:
        prefix, suffix = file_name.strip().split(".")
        sketch_type = ""
        if suffix == "fna":
            sketch_type = "dna"
        elif suffix == "faa":
            sketch_type = "protein"
        else:
            raise Exception(f"Unknown extension {suffix}.")
        print(f"Sketching the entries in the file {file_name}")
        make_sketches(ksize, scale_factor, file_name, sketch_type, out_dir, per_record=True)


if __name__ == "__main__":
    main()

