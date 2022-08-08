#!/usr/bin/env python
import argparse
import os
from Bio import SeqIO
import pandas as pd
import screed


def parse_args():
    parser = argparse.ArgumentParser(description="This script will take a simulation and a reference database of genes"
                                                 "and find the location of the genes in the simulation. The output will"
                                                 "include a bunch of helpful information about exactly what genes"
                                                 "were found in the simulation.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--database_dir", type=str, help="The location of the reference database")
    parser.add_argument("--simulation", type=str, help="The BBMB simulation to analyze")
    parser.add_argument("--output_file", type=str, help="The name of the output file")
    return parser.parse_args()


def main():
    # parse the arguments
    args = parse_args()
    database_dir = args.database_dir
    simulation_file = args.simulation
    output_file = args.output_file
    # check that the database directory exists and is not empty
    if not os.path.isdir(database_dir):
        print("Database directory is not a valid directory!")
        exit(-1)
    if len(os.listdir(database_dir)) == 0:
        print("Database directory is empty!")
        exit(-1)
    # check that the simulation file exists
    if not os.path.isfile(simulation_file):
        print("Simulation file is not a valid file!")
        exit(-1)
    # check that the output file doesn't already exist
    if os.path.isfile(output_file):
        print("Output file already exists! Overwriting...")
        os.remove(output_file)
    # get the list of genomes in the database
    genome_dir_names = [x[0] for x in os.walk(database_dir)][1:]
    genome_names = [os.path.split(x)[-1] for x in genome_dir_names]
    # find the mapping files in each of the genome dirs
    mapping_files = []
    for i in range(len(genome_dir_names)):
        mapping_file = os.path.join(genome_dir_names[i], f"{genome_names[i]}_mapping.csv")
        if not os.path.isfile(mapping_file):
            print(f"Mapping file {mapping_file} does not exist!")
            exit(-1)
        else:
            mapping_files.append(mapping_file)

    # import all the mapping files
    mapping_dfs = []
    for mapping_file in mapping_files:
        mapping_df = mapping_dfs.append(pd.read_csv(mapping_file))
    # merge all the mapping files into one dataframe
    mapping_df = pd.concat(mapping_dfs)

    # import the simulation file
    sequence_headers = []
    with screed.open(simulation_file) as seqfile:
        for read in seqfile:
            sequence_headers.append(read.name)

    # parse the headers
    # From Mahmudur, the header format is (separated by underscores):
    # SYN
    # id (0 indexed)
    # start (0 indexed)
    # end (0 indexed, the read includes the character at end position)
    # insert - not sure, but from the code, seems like this is the distance of the paired reads. This is meaningless when we use paired=false
    # strand
    # bbstart - from the code, it looks like this is some internal coordinate used by bbtools
    # bbchrom - again, from code, seems like internal representation of chromosomes inside of bbtools’ code
    simulation_read_locs = []
    for header in sequence_headers:
        header_split = header.strip().split("_")
        start = int(header_split[2])
        end = int(header_split[3])
        simulation_read_locs.append((start, end))
    # sanity check the read lengths
    #read_lengths = [x[1] - x[0] for x in simulation_read_locs]
    # make unique the read lengths
    #read_lengths = list(set(read_lengths))
    #print(read_lengths)




if __name__ == "__main__":
    main()



# ./find_genes_in_sim.py --database_dir ../test_data/input/reference_genomes/ --simulation ../test_data/output/simulated_metagenome.fq --output_file ../test_data/output/found_genes.csv