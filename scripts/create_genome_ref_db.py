#!/usr/bin/env python
import argparse
import os
import pandas as pd
from os import listdir
from os.path import isfile, join

def parse_args():
    parser = argparse.ArgumentParser(description="This script will take a directory as input where a list of genomes are "
                                                 " present. In each of the genome directory, there should be fna files "
                                                 "for each genome. These will be unioned into a single reference database, "
                                                 "from which metagenomes will be simulated. The number of genomes to be put in the"
                                                 " database is defined by the num_genomes arguments. The first that many genomes are used.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_dir", type=str, help="The full path to the directory that contains the genomes.")
    parser.add_argument("db_name", type=str, help="Name of the fasta database")
    parser.add_argument("num_genomes", type=int, help="Number of genomes to be put in the database")
    return parser.parse_args()


def main():
    args = parse_args()
    genome_path = args.in_dir
    fasta_filename = args.db_name
    num_genomes = args.num_genomes

    if not os.path.isdir(genome_path):
        print("Genome directory is not a valid directory!")
        exit(-1)

    # output fasta file
    fasta_file = open(fasta_filename, 'w')

    # iterate over the genomes
    genome_dir_names = [x[0] for x in os.walk(genome_path)][1:]
    if num_genomes > len(genome_dir_names):
        print("Given genome path: " + str(genome_path))
        print("Not enough genomes in the directory. Revise the num_genomes argument.")
        exit(-1)

    for genome_dir in genome_dir_names[:num_genomes]:
        genome_name = genome_dir.split('/')[-1]
        # find which of them is the genomic fna, but not cds or rna
        genome_dir_files = [f for f in listdir(genome_dir) if isfile(join(genome_dir, f))]
        reference_genome = None
        for file in genome_dir_files:
            if file.endswith("genomic.fna") and "rna_from_genomic" not in file and "cds_from_genomic" not in file:
                reference_genome = file
                reference_genome_path = os.path.join(genome_dir, reference_genome)
                break
        if reference_genome is None:
            print(f"Could not find reference genome for {genome_name}. Skipping...")
        else:
            # write the reference genome to the fasta file
            with open(reference_genome_path, 'r') as f:
                fasta_file.write(f.read())
            print(f"Wrote reference genome for {genome_name} to {fasta_filename}")
    fasta_file.close()
    print('DONE!')

if __name__ == "__main__":
    main()
