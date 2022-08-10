#!/usr/bin/env python
import argparse
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="This script will take a directory as input where a list of genomes are "
    " present. In each of the genome directory, there should be a mapping file present with the same name as the genome "
    " and with the extension csv. These csv files should be created using the script make_mapping_file.py. This script will then parse all of these csv "
    " files, and then create a reference database of all the genes present in those mapping files. The output of this script will be a fasta file, "
    " with one sequence for every gene.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_dir", type=str, help="The full path to the directory that contains the genomes.")
    parser.add_argument("db_name", type=str, help="Name of the fasta database")
    return parser.parse_args()

def main():
    args = parse_args()
    genome_path = args.in_dir
    fasta_filename = args.db_name

    if not os.path.isdir(genome_path):
        print("Genome directory is not a valid directory!")
        exit(-1)

    # keep record of all gene names to identify duplicates
    all_gene_names = set()

    # output fasta file
    fasta_file = open(fasta_filename, 'w')

    # iterate over the genomes
    genome_dir_names = [x[0] for x in os.walk(genome_path)][1:]
    for genome_dir in genome_dir_names:
        genome_name = genome_dir.split('/')[-1]
        mapping_filename = genome_name + "_mapping.csv"

        # check if mapping file exists
        if not os.path.exists(genome_dir + '/' + mapping_filename):
            print('Mapping file for genome ' + genome_name + ' was not found!')
            exit(-1)

        print('Handling genome: ' + genome_name)

        # identify repeating gene name (locus tag)
        df = pd.read_csv(genome_dir + '/' + mapping_filename)
        for gene_name in df['gene_name'].to_list():
            if gene_name in all_gene_names:
                print('ALERT: Gene with locus tag ' + gene_name + ' was found more than once!')
            else:
                all_gene_names.add(gene_name)

        # iterate over all recorded gene names, write in fasta
        for index, row in df.iterrows():
            genome_name = row['genome_name']
            assembly_id = row['assembly_id']
            gene_name = row['gene_name']
            protein_id = row['protein_id']
            contig_id = row['contig_id']
            start_position = row['start_position']
            end_position = row['end_position']
            strand = row['strand']
            aa_sequence = row['aa_sequence']

            header = '>' + gene_name + '|' + protein_id + '|' + assembly_id + '|' + contig_id + '|' + str(start_position) + '|' + str(end_position)
            fasta_file.write(header + '\n')
            fasta_file.write(aa_sequence + '\n')

    fasta_file.close()
    print('DONE!')

if __name__ == "__main__":
    main()
