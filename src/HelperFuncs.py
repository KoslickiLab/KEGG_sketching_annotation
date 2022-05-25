import os
import subprocess
import re
from collections import Counter
import numpy as np

bbtools_loc = os.path.abspath("../utils/bbmap")


def run_simulation(reference_file, out_file, num_reads, len_reads=150, noisy=False):
    """
    This function runs a simulation using bbtools "randomreads.sh"
    :param reference_file: The input sequences from which to make a metagenome TODO: make this auto downsample since currently it uses the whole set
    :param out_file: the name of the output simulation (must be a FASTQ file, so ending in fq or fastq)
    :param num_reads: number of reads to simulate
    :param len_reads: how long the reads are (default is 150bp)
    :param noisy: flag if you want noise injected to the simulation
    :return:
    """
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
        snprate = .01
        insrate = .01
        delrate = .01
        subrate = .01
        nrate = .01
        cmd += f"snprate={snprate} insrate={insrate} delrate={delrate} subrate={subrate} nrate={nrate} "
    cmd += f"ref={reference_file} out={out_file} reads={num_reads} length={len_reads} "
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    return


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
    print(f"I found {len(matches_tally)} unique matches totalling {np.sum(list(matches_tally.values()))} total matches "
          f"with the most frequent one occurring {np.max(list(matches_tally.values()))} times")
    return matches_tally


def make_sketches(ksize, scale_factor, file_name, sketch_type, out_dir, per_record=False):
    """
    This helper function will create the signature/sketches using sourmash
    :param ksize: the k-size to use
    :param scale_factor: the denominator of the scale factor to use (so >=1)
    :param file_name: the file to sketch
    :param sketch_type: amino acid (aa) or nucleotide (nt)
    :param out_dir: Where to write the signature
    :param per_record: If you want sketches of each entry in the fasta file, or just of the full fasta file (default: False)
    :return: None
    """
    out_file = f"{os.path.join(out_dir, os.path.basename(file_name))}_k_{ksize}_scale_{scale_factor}.sig"
    if per_record:
        cmd = f"sourmash sketch {sketch_type} -p k={ksize},scaled={scale_factor},abund -o {out_file} --singleton {file_name}"
    else:
        cmd = f"sourmash sketch {sketch_type} -p k={ksize},scaled={scale_factor},abund -o {out_file} --singleton {file_name}"
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)


# TODO: make a sourmash gather helper function
# TODO: make binary metric measures