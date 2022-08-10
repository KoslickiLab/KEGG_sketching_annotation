import os
import subprocess
import re
from collections import Counter
import numpy as np
import pandas as pd
import pathlib
import tempfile
from os.path import exists
import multiprocessing

bbtools_loc = os.path.abspath("../utils/bbmap")
diamond_loc = os.path.abspath("../utils/")


def run_simulation(reference_file, out_file, num_reads, len_reads=150, noisy=False, num_orgs=250):
    """
    This function runs a simulation using bbtools "randomreads.sh"
    :param reference_file: The input sequences from which to make a metagenome TODO: make this auto downsample since currently it uses the whole set
    :param out_file: the name of the output simulation (must be a FASTQ file, so ending in fq or fastq)
    :param num_reads: number of reads to simulate
    :param len_reads: how long the reads are (default is 150bp)
    :param noisy: flag if you want noise injected to the simulation
    :param num_orgs: specify the number of organisms/genes/proteins/etc. to include in the simulation
    :return: None
    """
    # First subsample the database so there are only num_orgs in the new reference file
    with tempfile.NamedTemporaryFile(suffix=pathlib.Path(reference_file).suffix) as subsample_ref_file:
    #with open('/tmp/test.faa', 'w') as subsample_ref_file:
        # do the subsampling
        cmd = f"{bbtools_loc}/./reformat.sh in={reference_file} out={subsample_ref_file.name} ow=t " \
              f"samplereadstarget={num_orgs} ignorejunk=t iupacToN=f crashjunk=f fixjunk=f"
        res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        if res.returncode != 0:
            raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
        # then generate the metagenome
        cmd = f"{bbtools_loc}/./randomreads.sh "
        # Things we always want set to true
        simple_names = "t"
        illumina_names = "f"
        overwrite = "t"
        metagenome = "t"
        banns = "t"
        # Note: reads by default are exactly 150bp long, not paired
        cmd += f"simplenames={simple_names} overwrite={overwrite} illuminanames={illumina_names} metagenome={metagenome} banns={banns} "
        if noisy:
            # TODO: hard code these for now, look up realistic values later
            snprate = 0.01
            insrate = 0.01
            delrate = 0.01
            subrate = 0.01
            nrate = 0.01
            cmd += f"snprate={snprate} insrate={insrate} delrate={delrate} subrate={subrate} nrate={nrate} "
        else:
            cmd += f"snprate=0 insrate=0 delrate=0 subrate=0 nrate=0 maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 maxns=0 adderrors=f "
        cmd += f"ref={subsample_ref_file.name} out={out_file} reads={num_reads} length={len_reads} "
        res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        if res.returncode != 0:
            raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    return


def compute_rel_abundance(simulation_fq_file):
    """
    This function is deprecated!!

    This helper function will find the true relative abundances of the sequences in the simulation
    :param simulation_fq_file: The bbmap simulated metagenome
    :return: a Counter object that contains the sequence identifiers and their counts in the simulation
    """
    # This assumes that the sequences are labeled as _._[three lower case letters}:{mixed case}_{integer}|
    #regex = r'\_\.\_([a-z]{3}:[a-zA-Z]\w+)'
    # Looks like the labels changed in bbmap v3.0.0, this is now _._[two lower case letters]_{mixed case}_{integer}|
    regex = r'\_\.\_([a-zA-Z]{2}_[a-zA-Z]\w+)'
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
    :param sketch_type: amino acid (aa, protein) or nucleotide (nt, dna)
    :param out_dir: Where to write the signature
    :param per_record: If you want sketches of each entry in the fasta file, or just of the full fasta file (default: False)
    :return: None
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = f"{os.path.join(out_dir, os.path.basename(file_name))}_k_{ksize}_scale_{scale_factor}.sig"
    if sketch_type == 'aa' or sketch_type == 'protein':
        sketch_type = "protein"
    elif sketch_type == 'nt' or sketch_type == 'dna':
        sketch_type = "dna"
    # Don't include this last line so we can work with translated DNA->Protein as well
    #else:
    #    raise Exception(f"Sketch_type must be one of 'aa' or 'nt. I was given {sketch_type}")
    if per_record:
        cmd = f"sourmash sketch {sketch_type} -f -p k={ksize},scaled={scale_factor},abund -o {out_file} --singleton {file_name}"
    else:
        cmd = f"sourmash sketch {sketch_type} -f -p k={ksize},scaled={scale_factor},abund -o {out_file} {file_name}"
    res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    if res.returncode != 0:
        raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    return


def run_sourmash_gather(query, database, out_file, sketch_type, num_results=None, threshold_bp=1000, quiet=True):
    """
    This is a simple wrapper for sourmash gather. It is hard coded to ignore abundances, estimate the
    ani and ci, as well as not perform the prefetch steps.
    :param query: file containing the query sketch/signature
    :param database: file containing the database sketch/signature
    :param out_file: the output csv file with the results
    :param sketch_type: aa (for amino acid) or nt (for nucleotide)
    :param num_results: int, if you only want the top N results
    :param threshold_bp: int, stop the algorithm once the overlap is below this many base pairs
    :return:
    """
    ignore_abundance = True
    estimate_ani_ci = True
    no_prefetch = False  # if set to true, this will disable the prefetch step. This will result in a much slower execution since prefetch does a pre-filtering step to select only database entries relevant to the sample
    if sketch_type not in ['aa', 'nt', 'protein', 'dna']:
        raise Exception(f"sketch type must be one of aa or protein (amino acid) or nt or dna (nucleotide). Provided value was {sketch_type}")
    cmd = f"sourmash gather -o {out_file} "
    if ignore_abundance:
        cmd += "--ignore-abundance "
    if estimate_ani_ci:
        cmd += "--estimate-ani-ci "
    if no_prefetch:
        cmd += "--no-prefetch "
    if sketch_type == 'aa' or sketch_type == 'protein':
        cmd += "--protein "
    elif sketch_type == 'nt' or sketch_type == 'dna':
        cmd += "--dna "
    if num_results:
        cmd += f"--num-results {num_results} "
    if threshold_bp:
        cmd += f"--threshold-bp {threshold_bp} "
    if quiet:
        cmd += "-q "
    cmd += f"{query} {database}"
    res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    if res.returncode != 0:
        raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    return


def return_unique_gather_hits(gather_out_file):
    """
    Takes a sourmash gather csv and returns the unique hits/identifiers in it
    :param gather_out_file: csv file from sourmash gather -o
    :return: a list of unique gene identifiers
    """
    if not os.path.exists(gather_out_file):
        raise Exception(f"File {gather_out_file} does not exist")
    df = pd.read_csv(gather_out_file)
    names = list(df['name'])
    name_ids = [x.split('|')[0] for x in names]
    name_ids_unique = set(name_ids)
    return list(name_ids_unique)


def calc_binary_stats_sourmash(simulation_file, gather_out_file):
    """
    This function is deprecated!!

    This function takes the simulation fastq file and the gather csv out file and calculates
    binary statistics from it: a dict with keys TP, FP, FN, precision, recall, F1.
    If you provide it a *.abund file, it just reads it as is
    :param simulation_file: Fastq file that contains the simulated reads, or relative abundances already
    :param gather_out_file: the results of running sourmash gather on those simulated reads
    :return: dict
    """
    # If the gt results are precomputed, just read them in
    simulation_gene_ids = set()
    ext = pathlib.Path(simulation_file).suffix
    # If the abund was passed, just read it in
    if ext == '.abund':
        with open(f"{simulation_file}", 'r') as fid:
            for line in fid.readlines():
                ident, count = line.strip().split('\t')
                simulation_gene_ids.add(ident)
    elif ext == '.fq':
        simulation_gene_ids = set(compute_rel_abundance(simulation_file).keys())
    else:
        raise Exception(f"Unknown file extension {ext}. Must be either fq or abund")
    gather_gene_ids = set(return_unique_gather_hits(gather_out_file))
    stats = dict()
    stats['TP'] = len(gather_gene_ids.intersection(simulation_gene_ids))
    stats['FP'] = len(gather_gene_ids.difference(simulation_gene_ids))
    stats['FN'] = len(simulation_gene_ids.difference(gather_gene_ids))
    stats['precision'] = stats['TP'] / float(stats['TP'] + stats['FP'])
    stats['recall'] = stats['TP'] / float(stats['TP'] + stats['FN'])
    if stats['TP']:
        stats['F1'] = 2 * stats['precision'] * stats['recall'] / float(stats['precision'] + stats['recall'])
    else:
        stats['F1'] = 0
    return stats


def check_extension(file_name):
    """
    Checks the file extension to see if it's protein or dna
    :param file_name: file name to check
    :return: 'protein' or 'dna'
    """
    suffix = pathlib.Path(file_name).suffix
    sketch_type = ""
    if suffix == ".fna":
        sketch_type = "dna"
    elif suffix == ".faa":
        sketch_type = "protein"
    else:
        raise Exception(f"Unknown extension {suffix}.")
    return sketch_type


def build_diamond_db(input_file, output_database):
    """
    This function is a simple wrapper for DIAMOND to create a reference database from protein sequences
    :param input_file: input reference FASTA file
    :param output_database: output database file (in DIAMOND binary format)
    :return: none
    """
    # DIAMOND will not automatically create folders
    if not exists(os.path.dirname(output_database)):
        os.makedirs(os.path.dirname(output_database))
    cmd = f"{diamond_loc}/./diamond makedb --in {input_file} -d {output_database}"
    res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    if res.returncode != 0:
        raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    return


def run_diamond_blastx(query_file, database_file, out_file, num_threads=multiprocessing.cpu_count()):
    """
    This is a simple wrapper to take a metagenome/query file `query_file` that contains DNA sequences,
    translate it to protein space, and then align against the database file.
    :param query_file: Input FASTA/Q query file
    :param database_file: The database built with build_diamond_db
    :param out_file: The output tsv file. Format is:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    :param num_threads: Number of threads to run (default=number of CPU cores on the machine you are using)
    :return: none
    """
    cmd = f"{diamond_loc}/./diamond blastx -d {database_file} -q {query_file} -o {out_file} -p {num_threads}"
    res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    if res.returncode != 0:
        raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    return


def parse_diamond_results(matches_file):
    """
    This parses the DIAMOND output csv file and returns some values about the results.
    :param matches_file: the output csv file from DIAMOND
    :return: 3-tuple: a) the set of gene identifiers that DIAMOND predicted to be in the sample
    b) the number of correct alignments (diamond aligned the read to the correct reference sequence)
    c) the number of incorrect alignments  (diamond aligned the read to the wrong reference sequence)
    """
    df = pd.read_csv(matches_file, sep='\t', header=0,
                       names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                              "send", "evalue", "bitscore"])
    query_names = list(df['qseqid'])
    regex = r'\_\.\_([a-z]{3}:[a-zA-Z]\w+)'
    regexc = re.compile(regex)
    query_ids = [re.findall(regexc, x)[0] for x in query_names]
    ref_ids = [x.split('|')[0] for x in df['sseqid']]
    if len(query_ids) != len(ref_ids):
        raise Exception(f"Something went wrong: there are {len(query_ids)} query ids, but {len(ref_ids)} ref ids.")
    inferred_ids = set()
    n_correct_alignments = 0
    n_incorrect_alignments = 0
    for true_id, inf_id in zip(query_ids, ref_ids):
        inferred_ids.add(inf_id)
        if true_id == inf_id:
            n_correct_alignments += 1
        else:
            n_incorrect_alignments += 1
    return inferred_ids, n_correct_alignments, n_incorrect_alignments


def calc_binary_stats_diamond(simulation_file, matches_file):
    """
    This function is deprecated!!

    This calculates the binary statistics (from a pure "gene present/absent" perspective) for the performance of DIAMOND
    on simulated data
    :param simulation_file: the simulation fastq file on which diamond was run
    :param matches_file: the output from diamond
    :return: dict (with stats in it)
    """
    # If the gt results are precomputed, just read them in
    simulation_gene_ids = set()
    ext = pathlib.Path(simulation_file).suffix
    # If the abund was passed, just read it in
    if ext == '.abund':
        with open(f"{simulation_file}", 'r') as fid:
            for line in fid.readlines():
                ident, count = line.strip().split('\t')
                simulation_gene_ids.add(ident)
    elif ext == '.fq':
        simulation_gene_ids = set(compute_rel_abundance(simulation_file).keys())
    else:
        raise Exception(f"Unknown file extension {ext}. Must be either fq or abund")
    diamond_gene_ids, n_correct_alignments, n_incorrect_alignments = parse_diamond_results(matches_file)
    diamond_gene_ids = set(diamond_gene_ids)
    stats = dict()
    stats['TP'] = len(diamond_gene_ids.intersection(simulation_gene_ids))
    stats['FP'] = len(diamond_gene_ids.difference(simulation_gene_ids))
    stats['FN'] = len(simulation_gene_ids.difference(diamond_gene_ids))
    stats['precision'] = stats['TP'] / float(stats['TP'] + stats['FP'])
    stats['recall'] = stats['TP'] / float(stats['TP'] + stats['FN'])
    if stats['TP']:
        stats['F1'] = 2 * stats['precision'] * stats['recall'] / float(stats['precision'] + stats['recall'])
    else:
        stats['F1'] = 0
    stats['Percent correct alignments'] = n_correct_alignments / float(n_correct_alignments + n_incorrect_alignments)
    stats['Total number of alignments'] = n_correct_alignments + n_incorrect_alignments
    cmd = f"wc -l {simulation_file}"
    res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    stats['Total number of sequences'] = int(res.stdout.split()[0]) / 4
    return stats


def parse_sourmash_results(gather_file):
    """
    This function will parse the output from sourmash gather and turn it into a functional profile that we can compare
    to the ground truth
    :param gather_file: the csv output from sourmash
    :return: a dataframe in the same format as that returned by find_genes_in_sim.py
    """
    pass



# TODO: calculate weighted stats. Need to understand what the difference columns in the sourmash gather results are actually returning