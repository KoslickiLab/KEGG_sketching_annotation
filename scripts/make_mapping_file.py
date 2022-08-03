#!/usr/bin/env python
import argparse
import os
from Bio import SeqIO
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="This script will create a mapping file for the genes present in\n"
                                                " a list of genomes. The genomes should be located in a directory.\n"
                                                " The directory is entered as an argument. The genomes should be in\n"
                                                " unzipped fna file format. The gbff files should be in the same directory.\n"
                                                " After finishing, the script will output one csv file per genome in that "
                                                " genome's directiry with the same name.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_dir", type=str, help="The full path to the directory that contains the genomes.")
    return parser.parse_args()

def main():
    # parse args
    args = parse_args()
    genome_path = args.in_dir

    if not os.path.isdir(genome_path):
        print("Not a valid directory!")
        exit(-1)

    genome_dir_names = [x[0] for x in os.walk(genome_path)][1:]
    for genome_dir in genome_dir_names:

        genome_name = genome_dir.split('/')[-1]
        fna_filename = None
        gbff_filename = None
        mapping_filename = genome_name + "_mapping.csv"

        for filename in list(os.walk(genome_dir))[0][2]:
            if filename.endswith("genomic.fna") and not "from" in filename:
                fna_filename = filename
            if filename.endswith("genomic.gbff"):
                gbff_filename = filename

        print('Handling genome: ' + genome_name)

        need_to_crash = False
        if fna_filename is None:
            print("Could not find fna file!")
            need_to_crash = True
        if gbff_filename is None:
            print("Could not find gbff file!")
            need_to_crash = True

        if need_to_crash:
            exit(-1)

        print('Working with fna file: ' + fna_filename + ' and gbff file: ' + gbff_filename)

        parsed_fna = SeqIO.parse(genome_dir + '/' + fna_filename, 'fasta')
        fasta_sequences = {}
        for rec in parsed_fna:
            fasta_sequences[ rec.id ] = rec.seq

        parsed_gbff = SeqIO.parse(genome_dir + '/' + gbff_filename, 'genbank')

        mapping_records = []
        for rec in parsed_gbff:
            for db in rec.dbxrefs:
                if db.startswith('Assembly'):
                    assembly_id = db.split(':')[1]
            for feature in rec.features:
                if feature.type == "CDS":
                    if 'protein_id' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
                        if 'translation' in feature.qualifiers:
                            if feature.location.start.position and feature.location.end.position:
                                gene_name = feature.qualifiers['locus_tag'][0]
                                start_position = feature.location.start.position
                                end_position = feature.location.end.position
                                strand = feature.strand  # Important to know if this sequence is read left to right or right to left
                                protein_id = feature.qualifiers['protein_id'][0]  # which protein this actually is
                                amino_acid = feature.qualifiers['translation'][0]  # the actual amino acid sequence that will be used to train the classifiers
                                mapping_records.append((genome_name, assembly_id, gene_name, protein_id, rec.id, start_position, end_position, strand, amino_acid, str(fasta_sequences[rec.id][start_position-1:end_position])))

        df = pd.DataFrame(mapping_records, columns=['genome_name', 'assembly_id', 'gene_name', 'protein_id', 'contig_id', 'start_position', 'end_position', 'strand', 'aa_sequence', 'nt_sequence'])
        df.to_csv(genome_dir + '/' + mapping_filename)

if __name__ == "__main__":
    main()
