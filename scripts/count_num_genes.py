#!/usr/bin/env python
import argparse
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="This script will take a directory as input where a list of genomes are "
    " present. In each of the genome directory, there should be a mapping file present with the same name as the genome "
    " and with the extension csv. These csv files should be created using the script make_mapping_file.py. This script will then "
    " count number of genes present in each genome.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_dir", type=str, help="The full path to the directory that contains the genomes.")
    parser.add_argument("num_genomes", type=int, help="Number of genomes to be put in the database")
    return parser.parse_args()

def main():
    args = parse_args()
    genome_path = args.in_dir
    num_genomes = args.num_genomes

    if not os.path.isdir(genome_path):
        print("Genome directory is not a valid directory!")
        exit(-1)

    # iterate over the genomes
    genome_dir_names = [x[0] for x in os.walk(genome_path)][1:]
    if num_genomes > len(genome_dir_names):
        print("Not enough genomes in the directory. Revise the num_genomes argument.")
        exit(-1)

    cumulative_count = 0
    for genome_dir in genome_dir_names[:num_genomes]:
        genome_name = genome_dir.split('/')[-1]
        mapping_filename = genome_name + "_mapping.csv"

        # check if mapping file exists
        if not os.path.exists(genome_dir + '/' + mapping_filename):
            print('Mapping file for genome ' + genome_name + ' was not found!')
            exit(-1)

        # identify repeating gene name (locus tag)
        df = pd.read_csv(genome_dir + '/' + mapping_filename)
        count = len(df['gene_name'].to_list())
        cumulative_count += count
        print(genome_name, str(count), str(cumulative_count) )

if __name__ == "__main__":
    main()
