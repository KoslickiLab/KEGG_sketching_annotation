#!/usr/bin/env python
import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess
import re
import numpy as np
from collections import Counter

bbtools_loc = os.path.abspath("../utils/bbmap")


def run_simulation(reference_file, out_file, num_reads, len_reads, noisy=False):
    cmd = f"{bbtools_loc}/./randomreads.sh "
    # Things we always want set to true
    simple_names = "t"
    illumina_names = "f"
    overwrite = "t"
    metagenome = "t"
    # Note: reads by default are exactly 150bp long, not paired
    cmd += f"simplenames={simple_names} overwrite={overwrite} illuminanames={illumina_names} metagenome={metagenome} "
    #cmd += f"metagenome={metagenome} "
    if noisy:
        # TODO: hard code these for now, look up realistic values later
        snprate = .51
        insrate = .01
        delrate = .01
        subrate = .01
        nrate = .01
        cmd += f"snprate={snprate} insrate={insrate} delrate={delrate} subrate={subrate} nrate={nrate} "
    cmd += f"ref={reference_file} out={out_file} reads={num_reads} length={len_reads} "
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)


def compute_rel_abundance(simulation_fq_file):
    """
    This helper function will find the true relative abundances of the sequences in the simulation
    :param simulation_fq_file: The bbmap simulated metagenome
    :return: a Counter object that contains the sequence identifiers and their counts in the simulation
    """
    # This assumes that the sequences are labeled as _._[three lower case letters}:{mixed case}_{integer}|
    regex = r'\_\.\_([a-z]{3}:[a-zA-Z]\w+)'
    contents = open(simulation_fq_file, 'r').read()
    matches = re.findall(regex, contents)
    matches_tally = Counter(matches)
    print(f"I found {len(matches_tally)} unique matches totalling {np.sum(list(matches_tally.values()))} total matches")


def main():
    parser = argparse.ArgumentParser(description="This script uses BBTools randomreads.sh to simulate reads from a"
                                                 "metagenome.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reference_file', type=str, help="Database that you want to generate reads from.")
    parser.add_argument('-o', '--out_file', type=str, help="The output simulated metagenome.")
    parser.add_argument('-n', '--num_reads', type=int, help="The number of reads to simulate.")
    parser.add_argument('-l', '--len_reads', type=int, help="The length of the reads to generate.")
    parser.add_argument("--noisy", action='store_true', help="The full path to the directory that the files will be written")
    # parse the args
    args = parser.parse_args()
    reference_file = args.reference_file
    out_file = args.out_file
    ext = out_file.split(".")[-1]
    if ext not in ["fq", "fastq"]:
        raise Exception(f"Output file extension must be one of fq or fastq. Yours was {ext}")
    num_reads = args.num_reads
    len_reads = args.len_reads
    noisy = args.noisy
    # run the simulation
    run_simulation(reference_file, out_file, num_reads, len_reads, noisy=noisy)
    # For some really odd reason, some (but not all) of the underscores in the names are converted to left brackets {
    # So let's replace them back
    subprocess.run("sed -i 's/{/_/g' " + out_file, stdout=subprocess.PIPE, shell=True)



if __name__ == "__main__":
    main()

