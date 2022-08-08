#!/usr/bin/env python
import argparse
import os
from Bio import SeqIO
import pandas as pd


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
    genome_names = [x.split('/')[-1] for x in genome_dir_names]
