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

# TODO: this is very slow O(n^2), could make it faster by using dicts or the like


def interval_length(intvals):
    length = 0
    for intval in intvals:
        length += intval[1] - intval[0] + 1
    return int(length)


def create_ground_truth(found_genes_df):
    """
    This function will take the dataframe that contains all the genes found in the simulation and output
    a ground truth functional profile
    :param found_genes_df:
    :return: dataframe with columns: gene_name, nucleotide_overlap, median_coverage, mean_coverage, reads_mapped
    """
    unique_genes = found_genes_df.gene_name.unique()
    ground_truth_df = pd.DataFrame(columns=["gene_name", "gene_length", "nucleotide_overlap", "median_coverage", "mean_coverage", "reads_mapped"])
    ground_truth_df_data = {"gene_name": [], "gene_length": [], "nucleotide_overlap": [], "median_coverage": [], "mean_coverage": [], "reads_mapped": [], "num_nts_in_reads": []}
    # TODO: the iterations of this are slow. The iterations that have a large len(intervals) are particularly slow.
    # TODO: it might be the unioning of all these intervals that's slow, or the population of the coverage array
    # TODO: which is O(gene_length * len(intervals))
    itr = 0
    for gene in unique_genes:
        itr += 1
        #print(f"Analyzing gene: {itr}/{len(unique_genes)}")
        gene_df = found_genes_df[found_genes_df.gene_name == gene]
        # get overlap intervals
        interval_tuples = list(zip(gene_df.overlap_start, gene_df.overlap_end))
        # Convert these to intervals
        intervals = [interval[start, end] for start, end in interval_tuples]
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
        # FIXME: this is slow, it's O(gene_length * len(intervals))
        #for i in range(gene_length):
        #    coverage_array[i] = sum((i+gene_start) in intval for intval in intervals)  # sum the coverage for each position in the gene (the array is True or False if that i position is covered)
        for intval in intervals:
            # Some of the intervals are singletons, so check these
            left_end = int(intval.extrema[0][0])
            try:
                right_end = int(intval.extrema[1][0])
            except IndexError:
                right_end = left_end
            coverage_array[(left_end - gene_start):(right_end - gene_start)] += 1

        median_coverage = np.median(coverage_array)
        mean_coverage = np.mean(coverage_array)

        num_nts_in_reads = sum(gene_df.num_bases_overlap)

        # store the information in the dataframe
        ground_truth_df_data["gene_name"].append(gene)
        ground_truth_df_data["gene_length"].append(gene_length)
        ground_truth_df_data["nucleotide_overlap"].append(nucleotide_overlap)
        ground_truth_df_data["median_coverage"].append(median_coverage)
        ground_truth_df_data["mean_coverage"].append(mean_coverage)
        ground_truth_df_data["reads_mapped"].append(reads_mapped)
        ground_truth_df_data["num_nts_in_reads"].append(num_nts_in_reads)
    ground_truth_df = pd.DataFrame(ground_truth_df_data)
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
    elif start > end_pos:  # no overlap
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
    parser.add_argument("--num_genomes", type=int, help="Number of genomes to consider")
    return parser.parse_args()


def main():
    # parse the arguments
    args = parse_args()
    database_dir = args.database_dir
    simulation_file = args.simulation
    output_file = args.output_file
    num_genomes = args.num_genomes

    # check that the database directory exists and is not empty
    if not os.path.isdir(database_dir):
        print("Database directory is not a valid directory!")
        exit(-1)
    if len(os.listdir(database_dir)) == 0:
        print("Database directory is empty!")
        exit(-1)
    # check that there are enough genomes in the directory
    genome_dir_names = [x[0] for x in os.walk(database_dir)][1:]
    if num_genomes > len(genome_dir_names):
        print("Not enough genomes in the directory. Revise the num_genomes argument.")
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
    genome_dir_names = genome_dir_names[:num_genomes]
    genome_names = [os.path.split(x)[-1] for x in genome_dir_names]
    # find the mapping files in each of the genome dirs
    print("Finding mapping files...")
    mapping_files = []
    for i in range(len(genome_dir_names)):
        mapping_file = os.path.join(genome_dir_names[i], f"{genome_names[i]}_mapping.csv")
        if not os.path.isfile(mapping_file):
            print(f"Mapping file {mapping_file} does not exist!")
            exit(-1)
        else:
            mapping_files.append(mapping_file)

    # import all the mapping files
    print("Importing mapping files...")
    mapping_dfs = []
    for mapping_file in mapping_files:
        mapping_df = mapping_dfs.append(pd.read_csv(mapping_file))
    # merge all the mapping files into one dataframe
    mapping_df = pd.concat(mapping_dfs)
    print("num of records in all mapping files: " + str(len(mapping_df.index)))
    print(mapping_df.sample(5))
    #print(mapping_df)

    # import the simulation file
    print("Importing simulation file...")
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
    print("Parsing headers...")
    simulation_df = pd.DataFrame(columns=["contig_id", "start", "end"])
    simulation_data = {"contig_id": [], "start": [], "end": []}
    for header in sequence_headers:
        header_split = header.split("_")
        start = int(header_split[2])
        end = int(header_split[3])
        contig_id = "_".join(header_split[9:11])  # this assumes that the format of the simulation will not change
        if "_".join(header_split[9:11])[0] == '$':
            contig_id = contig_id.split("$")[1]
        else:
            contig_id = contig_id.split("$")[0]
        # add the row to the dataframe
        simulation_data["contig_id"].append(contig_id)
        simulation_data["start"].append(start)
        simulation_data["end"].append(end)
    simulation_df = pd.DataFrame(simulation_data)
    print("num of records in simulation dataframe: " + str(len(simulation_df.index)))
    print(simulation_df.sample(5))


    # create dataframe for output
    output_df = pd.DataFrame(columns=["contig_id", "gene_name", "num_bases_overlap", "overlap_start", "overlap_end", "gene_start", "gene_end"])
    # iterate through the simulation dataframe

    # Create a couple of data structures that will make it faster to determine if a read overlaps a gene
    print("Creating data structures for overlaps...")
    contig_intervals = {}
    contig_intervals_union = {}
    for index, row in mapping_df.iterrows():
        contig_id = row["contig_id"]
        if contig_id not in contig_intervals:
            contig_intervals[contig_id] = {}
            contig_intervals_union[contig_id] = interval[0, 0]
        gene_name = row["gene_name"]
        if gene_name not in contig_intervals[contig_id]:
            contig_intervals[contig_id][gene_name] = (row["start_position"], row["end_position"])
            contig_intervals_union[contig_id] = contig_intervals_union[contig_id] | interval[row["start_position"], row["end_position"]]
        else:
            raise Exception("Duplicate gene name in mapping file!")

    print("Finding overlaps...")
    output_df_data = {"contig_id": [], "gene_name": [], "num_bases_overlap": [], "overlap_start": [], "overlap_end": [], "gene_start": [], "gene_end": []}
    # iterate through the all the reads in the simulation file looking for overlaps to genes
    for i in range(len(simulation_df)):
        #if i % 1000 == 0:
        #    print(f"On step: {i}/{len(simulation_df)}")
        # get the contig id
        contig_id = simulation_df.iloc[i]["contig_id"]
        # get the start and end positions of the read
        start = simulation_df.iloc[i]["start"]
        end = simulation_df.iloc[i]["end"]
        # check if the start or then end of the read is within a gene region. Skip it if not.
        if contig_id in contig_intervals:  # if there are no genes on this contig, skip it
            # This assumes that the reads will always be shorter than the genes
            if start not in contig_intervals_union[contig_id] and end not in contig_intervals_union[contig_id]:
                continue
        else:
            continue
        gene_names = contig_intervals[contig_id].keys()
        # iterate through the gene names
        if gene_names:
            for gene_name in gene_names:
                # get the start and end positions of the gene
                gene_start = contig_intervals[contig_id][gene_name][0]
                gene_end = contig_intervals[contig_id][gene_name][1]
                overlap_interval = interval_overlap((start, end), (gene_start, gene_end))
                if overlap_interval:
                    # get the protein id
                    # get the number of bases that overlap
                    num_bases_overlap = overlap_interval[1] - overlap_interval[0] + 1
                    # get the overlap start and end positions
                    overlap_start = overlap_interval[0]
                    overlap_end = overlap_interval[1]
                    # add the information to the output dataframe
                    output_df_data["contig_id"].append(contig_id)
                    output_df_data["gene_name"].append(gene_name)
                    output_df_data["num_bases_overlap"].append(num_bases_overlap)
                    output_df_data["overlap_start"].append(overlap_start)
                    output_df_data["overlap_end"].append(overlap_end)
                    output_df_data["gene_start"].append(gene_start)
                    output_df_data["gene_end"].append(gene_end)
    print("Putting data in dataframe...")
    output_df = pd.DataFrame(output_df_data)
    #output_df.to_csv('/home/dkoslicki/Documents/KEGG_sketching_annotation/scripts/temp.csv', index=False)

    # write the output file
    # output_df.to_csv(output_file, index=False)
    print("Summarizing gene coverage information...")
    ground_truth_df = create_ground_truth(output_df)
    print("Writing output file...")
    ground_truth_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()



# ./find_genes_in_sim.py --database_dir ../test_data/input/reference_genomes/ --simulation ../test_data/output/simulated_metagenome.fq --output_file ../test_data/output/found_genes.csv
# database_dir="test_data/input/reference_genomes"
# simulation_file="test_data/output/simulated_metagenome.fq"
