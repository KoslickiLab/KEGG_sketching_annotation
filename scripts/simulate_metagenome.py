#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
import platform
# for relative imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from src.HelperFuncs import run_simulation

bbtools_loc = os.path.abspath("../utils/bbmap")


def main():
    parser = argparse.ArgumentParser(description="This script uses BBTools randomreads.sh to simulate reads from a"
                                                 "metagenome.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference_file', type=str, help="Database that you want to generate reads from.")
    parser.add_argument('-o', '--out_file', type=str, help="The output simulated metagenome.")
    parser.add_argument('-n', '--num_reads', type=int, help="The number of reads to simulate.")
    parser.add_argument('-l', '--len_reads', type=int, help="The length of the reads to generate.", default=150)
    parser.add_argument('--num_orgs', type=int, help="The number of organisms/genes/proteins to include in the simulation.", default=250)
    parser.add_argument('--seed', type=int, help="The seed to run the simulation.", default=0)
    parser.add_argument("--noisy", action='store_true', help="If you want to inject noise in the simulated reads")

    # parse the args
    args = parser.parse_args()
    num_orgs = args.num_orgs
    reference_file = os.path.abspath(args.reference_file)
    out_file = os.path.abspath(args.out_file)
    ext = out_file.split(".")[-1]
    if ext not in ["fq", "fastq"]:
        raise Exception(f"Output file extension must be one of fq or fastq. Yours was {ext}")
    num_reads = args.num_reads
    if not num_reads:
        raise Exception("Must specify the number of reads via --num_reads")
    len_reads = args.len_reads
    if not len_reads:
        raise Exception("Must specify the length of the reads via --len_reads")
    noisy = args.noisy
    seed = args.seed

    # run the simulation
    run_simulation(reference_file, out_file, num_reads, len_reads, noisy=noisy, num_orgs=num_orgs, seed=seed)
    # For some really odd reason, some (but not all) of the underscores in the names are converted to left brackets {
    # So let's replace them back
    if platform.system() == 'Darwin':
        subprocess.run("sed -i.bu 's/{/_/g' " + out_file, stdout=subprocess.PIPE, shell=True)
    else:
        subprocess.run("sed -i 's/{/_/g' " + out_file, stdout=subprocess.PIPE, shell=True)
    # This is now depreciated in favor of find_genes_in_sim.py
    #matches_tally = compute_rel_abundance(out_file)
    # Save this for later use
    # Sort it first, descending order
    #matches_tally_sorted = {k: v for k, v in sorted(matches_tally.items(), key=lambda item: -item[1])}
    #with open(f"{out_file}.abund", 'w') as fid:
    #    for key, value in matches_tally_sorted.items():
    #        fid.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    # Example of running this script
    # ./simulate_metagenome.py -r ../test_data/input/kegg_genes_KO.fna -o ../test_data/output/test_simulation.fq -n 1000 -l 150 --num_orgs 10
    main()
