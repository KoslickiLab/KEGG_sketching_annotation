#!/usr/bin/env python
import argparse
import os, sys
from Bio import SeqIO
import pandas as pd
import screed
import numpy as np
# add parent folder to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from utils.pyinterval.interval import interval


def interval_length(intvals):
    length = 0
    for intval in intvals:
        length += intval[1] - intval[0]
    return length


def create_ground_truth(found_genes_df):
    """
    This function will take the dataframe that contains all the genes found in the simulation and output
    a ground truth functional profile
    :param found_genes_df:
    :return: dataframe with columns: gene_name, nucleotide_overlap, median_coverage, mean_coverage, reads_mapped
    """
    unique_genes = found_genes_df.gene_name.unique()
    ground_truth_df = pd.DataFrame(columns=["gene_name", "nucleotide_overlap", "median_coverage", "mean_coverage", "reads_mapped"])
    for gene in unique_genes:
        gene_df = found_genes_df[found_genes_df.gene_name == gene]
        # get overlap intervals
        interval_tuples = list(zip(gene_df.overlap_start, gene_df.overlap_end))
        # Convert these to intervals
        intervals = [interval[start, end] for start, end in interval_tuples]
        #print(intervals)
        # take the union of all the intervals
        intervals_union = intervals[0]
        for intval in intervals[1:]:
            intervals_union = intervals_union | intval
        # this is the total number of gene nucleotides that have been covered by reads from the simulations
        nucleotide_overlap = interval_length(intervals_union)
        # Since I iterate over the reads, the number of entries in this data frame is the number of reads that has mapped to this gene
        reads_mapped = len(gene_df)
        # since this df is the subset corresponding to a single gene, I can pull off the length from the first one
        gene_start = gene_df.iloc[0].gene_start
        gene_end = gene_df.iloc[0].gene_end
        gene_length = gene_end - gene_start + 1
        # for each location in the gene, get the coverage information
        coverage_array = np.zeros(gene_length)
        for i in range(gene_length):
            coverage_array[i] = sum((i+gene_start) in intval for intval in intervals)
        median_coverage = np.median(coverage_array)
        mean_coverage = np.mean(coverage_array)
        # store the information in the dataframe
        ground_truth_df.loc[len(ground_truth_df)] = [gene, int(nucleotide_overlap), median_coverage, mean_coverage, int(reads_mapped)]
    return ground_truth_df


def interval_overlap(interval1, interval2):
    """
    This function will determine if two intervals overlap, and by how much.
    :param interval1: The first interval
    :param interval2: The second interval
    :return: True if the intervals overlap, False otherwise
    """
    start, end = interval1[0], interval1[1]
    start_pos, end_pos = interval2[0], interval2[1]
    if end < start_pos:  # no overlap
        return False
    elif end >= start_pos > start:  # right side overlap
        return start_pos, end
    elif start_pos <= end <= end_pos and start_pos <= start <= end_pos:  # proper subset overlap
        return start, end
    elif end > end_pos >= start >= start_pos:  # left side overlap
        return start, end_pos
    elif start > end_pos: # no overlap
        return False
    elif start < start_pos and end > end_pos:  # proper superset overlap
        return start_pos, end_pos
    else:
        print("Error in interval overlap")
        exit(-1)


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
    #print(mapping_df)

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
    simulation_df = pd.DataFrame(columns=["contig_id", "start", "end"])
    for header in sequence_headers:
        header_split = header.split("_")
        start = int(header_split[2])
        end = int(header_split[3])
        contig_id = "_".join(header_split[9:11])  # this assumes that the format of the simulation will not change
        contig_id = contig_id.split("$")[0]
        # add the row to the dataframe
        simulation_df.loc[len(simulation_df)] = [contig_id, start, end]

    # create dataframe for output
    output_df = pd.DataFrame(columns=["contig_id", "gene_name", "protein_id", "num_bases_overlap", "overlap_start", "overlap_end", "gene_start", "gene_end"])
    # iterate through the simulation dataframe
    for i in range(len(simulation_df)):
        # get the contig id
        contig_id = simulation_df.iloc[i]["contig_id"]
        # get the start and end positions of the read
        start = simulation_df.iloc[i]["start"]
        end = simulation_df.iloc[i]["end"]
        #print(f"Analyzing contig {contig_id}: start {start} end {end}")
        # get the mapping dataframe for the contig
        contig_mapping_df = mapping_df[mapping_df["contig_id"] == contig_id]
        # iterate through the mapping dataframe
        for j in range(len(contig_mapping_df)):
            # check if the read overlaps with this gene
            gene_start = contig_mapping_df.iloc[j]["start_position"]
            gene_end = contig_mapping_df.iloc[j]["end_position"]
            overlap_interval = interval_overlap((start, end), (gene_start, gene_end))
            if overlap_interval:
                # get the gene name
                gene_name = contig_mapping_df.iloc[j]["gene_name"]
                # get the protein id
                protein_id = contig_mapping_df.iloc[j]["protein_id"]
                # get the number of bases that overlap
                num_bases_overlap = overlap_interval[1] - overlap_interval[0] + 1
                # get the overlap start and end positions
                overlap_start = overlap_interval[0]
                overlap_end = overlap_interval[1]
                # add the row to the dataframe
                output_df.loc[len(output_df)] = [contig_id, gene_name, protein_id, num_bases_overlap, overlap_start, overlap_end, gene_start, gene_end]
                #print(f"Gene start {gene_start} end {gene_end}, read start {start} end {end}")
                #print(f"Overlap start {overlap_start} end {overlap_end}")
    # write the output file
    # output_df.to_csv(output_file, index=False)
    ground_truth_df = create_ground_truth(output_df)
    ground_truth_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()



# ./find_genes_in_sim.py --database_dir ../test_data/input/reference_genomes/ --simulation ../test_data/output/simulated_metagenome.fq --output_file ../test_data/output/found_genes.csv
# database_dir="test_data/input/reference_genomes"
# simulation_file="test_data/output/simulated_metagenome.fq"