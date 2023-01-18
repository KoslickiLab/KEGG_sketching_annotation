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
import matplotlib.pyplot as plt
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
bbtools_loc = os.path.abspath(f"{THIS_DIR}/../utils/bbmap")
diamond_loc = os.path.abspath(f"{THIS_DIR}/../utils/")


def run_simulation(reference_file, out_file, num_reads, len_reads=150, noisy=False, num_orgs=250, seed=0):
    """
    This function runs a simulation using bbtools "randomreads.sh"
    :param reference_file: The input sequences from which to make a metagenome
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
            snprate = 0.001
            insrate = 0.001
            delrate = 0.001
            subrate = 0.001
            nrate = 0.001
            cmd += f"snprate={snprate} insrate={insrate} delrate={delrate} subrate={subrate} nrate={nrate} "
        else:
            cmd += f"snprate=0 insrate=0 delrate=0 subrate=0 nrate=0 maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 maxns=0 adderrors=f "
        cmd += f"ref={subsample_ref_file.name} out={out_file} reads={num_reads} length={len_reads} seed={seed}"
        res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        if res.returncode != 0:
            raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    return


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
    ignore_abundance = False
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
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore.
    Or put more verbosely: # Fields: Query ID, Subject ID, Percentage of identical matches, Alignment length,
    Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query,
    Start of alignment in subject, End of alignment in subject, Expected value, Bit score
    :param num_threads: Number of threads to run (default=number of CPU cores on the machine you are using)
    :return: none
    """
    cmd = f"{diamond_loc}/./diamond blastx -d {database_file} -q {query_file} -o {out_file} -p {num_threads}"
    #res = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    res = subprocess.run(cmd.split(' '))
    if res.returncode != 0:
        print(cmd)
        raise Exception(f"The command {cmd} exited with nonzero exit code {res.returncode}")
    # put the header into the file
    header = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
    with open(out_file, 'r') as original:
        data = original.read()
    with open(out_file, 'w') as modified:
        modified.write(f"{header}\n" + data)
    return


def parse_diamond_results(matches_file):
    """
    This parses the DIAMOND output csv file and returns some values about the results.
    :param matches_file: the output csv file from DIAMOND
    :return: 3-tuple: a) the set of gene identifiers that DIAMOND predicted to be in the sample
    b) the number of correct alignments (diamond aligned the read to the correct reference sequence)
    c) the number of incorrect alignments  (diamond aligned the read to the wrong reference sequence)
    """
    df = pd.read_csv(matches_file, sep='\t')
    ref_ids = [x.split('|')[0] for x in df['sseqid']]
    ref_ids_tally = Counter(ref_ids)
    # also get the gene lengths
    gene_lengths = [int(x.split('|')[-1]) - int(x.split('|')[-2]) + 1 for x in df['sseqid']]
    bit_scores = df['bitscore']
    ref_id_to_bit_score = dict()
    # add all the bit scores up for each reference id
    for i in range(len(ref_ids)):
        if ref_ids[i] not in ref_id_to_bit_score:
            ref_id_to_bit_score[ref_ids[i]] = bit_scores[i]
        else:
            ref_id_to_bit_score[ref_ids[i]] += bit_scores[i]
    # divide by the total number of reads to get the average bit score for each reference id
    for ref_id in ref_id_to_bit_score:
        ref_id_to_bit_score[ref_id] /= ref_ids_tally[ref_id]
    id_to_length = dict(zip(ref_ids, gene_lengths))
    data = {"name": [], "num_reads": [], "ave_bit_score": [], "gene_length": []}
    for ref_id in ref_ids_tally:
        data["name"].append(ref_id)
        data["num_reads"].append(ref_ids_tally[ref_id])
        data["gene_length"].append(id_to_length[ref_id])
        data["ave_bit_score"].append(ref_id_to_bit_score[ref_id])
    return pd.DataFrame(data)


def calculate_diamond_performance(diamond_file, ground_truth_file, filter_threshold=0, bitscore_threshold=0):
    """
    This function will parse the output from Diamond and turn it into a functional profile, and then calculate
    performance statistics similar to calculate_sourmash_performance
    :param diamond_file: the matches.csv file that's output from ./classify_diamond.py
    :param ground_truth_file: The ground_truth.csv file that's output from find_genes_in_sim.py
    :param filter_threshold: Ignore genes in the ground truth that have fewer than filter_threshold reads mapping
    to them. Also remove these genes from Diamond.
    :param bitscore_threshold: Ignore Diamond entries that have a bit score less than bitscore_threshold
    :return: dataframe
    """
    # check if files exist
    if not exists(diamond_file):
        raise Exception(f"{diamond_file} does not exist")
    if not exists(ground_truth_file):
        raise Exception(f"{ground_truth_file} does not exist")
    # parse the DIAMOND output file
    ddf = parse_diamond_results(diamond_file)
    # parse the ground truth file
    gdf = pd.read_csv(ground_truth_file)
    # filter out genes that are too short
    genes_removed = set(gdf[gdf['reads_mapped'] < filter_threshold].gene_name)
    gdf = gdf[~gdf['gene_name'].isin(genes_removed)]
    ddf = ddf[~ddf['name'].isin(genes_removed)]
    # filter out genes that have a bit score less than the threshold
    ddf = ddf[ddf['ave_bit_score'] >= bitscore_threshold]
    # Get gene names
    ggene_names = gdf.gene_name
    dgene_names = ddf.name
    # sort everything by gene name
    gdf = gdf.sort_values(by='gene_name')
    # get binary data frames
    ddf_TP = ddf[ddf['name'].isin(ggene_names)]  # true positives in sourmash
    ddf_TP = ddf_TP.sort_values(by='name')  # sort by gene name
    ddf_FP = ddf[~ddf['name'].isin(ggene_names)]  # false positives in sourmash
    ddf_FN = gdf[~gdf['gene_name'].isin(dgene_names)]  # false negatives (ones in ground truth that aren't in sourmash)
    gdf_TP = gdf[gdf['gene_name'].isin(dgene_names)]  # subset the ground truth to concentrate on the ones in the gather results
    gdf_TP = gdf_TP.sort_values(by='gene_name')
    metrics = ['TP', 'FP', 'FN', 'precision', 'recall', 'F1', 'corr_reads_mapped', 'L1_reads_mapped',
               'corr_reads_mapped_div_gene_len', 'L1_reads_mapped_div_gene_len', 'percent_correct_predictions', 'total_number_of_predictions']
    # calculate the performance metrics
    performance = dict()
    performance['TP'] = len(ddf_TP)
    performance['FP'] = len(ddf_FP)
    performance['FN'] = len(ddf_FN)
    if float(performance['TP'] + performance['FP']) > 0:
        performance['precision'] = performance['TP'] / float(performance['TP'] + performance['FP'])
    else:
        performance['precision'] = 0
    if float(performance['TP'] + performance['FN']) > 0:
        performance['recall'] = performance['TP'] / float(performance['TP'] + performance['FN'])
    else:
        performance['recall'] = 0
    if performance['TP']:
        performance['F1'] = 2 * performance['precision'] * performance['recall'] / float(
            performance['precision'] + performance['recall'])
    else:
        performance['F1'] = 0
    performance['corr_reads_mapped'] = np.corrcoef(ddf_TP['num_reads'], gdf_TP['reads_mapped'])[0][1]
    ddf_TP_vec = np.array(ddf_TP['num_reads'].values)
    ddf_TP_vec = ddf_TP_vec / np.sum(ddf_TP_vec)
    gdf_TP_vec = np.array(gdf_TP['reads_mapped'].values)
    gdf_TP_vec = gdf_TP_vec / np.sum(gdf_TP_vec)
    performance['L1_reads_mapped'] = np.sum(np.abs(ddf_TP_vec - gdf_TP_vec))
    reads_mapped_div_gene_len = np.array(gdf_TP['reads_mapped'] / gdf_TP['gene_length'])
    reads_mapped_div_gene_len = reads_mapped_div_gene_len / np.sum(reads_mapped_div_gene_len)
    pred_reads_mapped_div_gene_len = np.array(ddf_TP['num_reads'] / ddf_TP['gene_length'])
    pred_reads_mapped_div_gene_len = pred_reads_mapped_div_gene_len / np.sum(pred_reads_mapped_div_gene_len)
    performance['corr_reads_mapped_div_gene_len'] = np.corrcoef(pred_reads_mapped_div_gene_len, reads_mapped_div_gene_len)[0][1]
    L1_average_abund_reads_mapped_div_gene_len = np.sum(np.abs(pred_reads_mapped_div_gene_len - reads_mapped_div_gene_len))
    performance['L1_reads_mapped_div_gene_len'] = L1_average_abund_reads_mapped_div_gene_len
    if float(len(ddf_TP) + len(ddf_FP)) > 0:
        performance['percent_correct_predictions'] = len(ddf_TP) / float(len(ddf_TP) + len(ddf_FP))
    else:
        performance['percent_correct_predictions'] = 0
    performance['total_number_of_predictions'] = len(ddf_TP) + len(ddf_FP)
    # also record what filter threshold was used
    performance['filter_threshold'] = filter_threshold
    # put this in a dataframe
    performance_df = pd.DataFrame(performance, index=[0])
    return performance_df

def calculate_sourmash_performance(gather_file, ground_truth_file, filter_threshold):
    """
    This function will parse the output from sourmash gather and turn it into a functional profile that we can compare
    to the ground truth.
    From the check_sourmash_correlation method, it appears that:
     f_unique_weighted correlates with reads mapped, median, mean coverage, and nucleotide_overlap
     reads mapped / gene length correlates with sourmash's average_abund and median_abund: corr=0.9957443813164701

    :param gather_file: the csv output from sourmash
    :param ground_truth_file: the ground truth file (output from find_genes_in_sim.py)
    :param filter_threshold: ignore ground truth genes that have less than this many reads in the ground truth
    :return: a dataframe containing the performance characteristics of the sourmash results
    """
    if not os.path.exists(gather_file):
        raise Exception(f"File {gather_file} does not exist")
    if not os.path.exists(ground_truth_file):
        raise Exception(f"File {ground_truth_file} does not exist")

    sourmash_rel_abund_col = 'f_unique_weighted'
    ground_truth_rel_abund_col = 'reads_mapped'

    # s prefix is for "sourmash" while g prefix is for "ground truth"
    sdf = pd.read_csv(gather_file)
    gdf = pd.read_csv(ground_truth_file)
    # sort the ground truth by gene name, this will be the order that we stick with
    gdf = gdf.sort_values(by='gene_name')

    # grab the sourmash infered gene names
    sgene_names = list(sdf['name'])
    sgene_names = [x.split('|')[0] for x in sgene_names]
    # replace the name column with the sourmash gene names
    sdf['name'] = sgene_names
    ggene_names = list(gdf['gene_name'])
    greads_mapped = np.sum(list(gdf['reads_mapped']))

    # remove the infrequent genes from the ground truth and from sourmash
    genes_removed = set(gdf[gdf['reads_mapped'] < filter_threshold].gene_name)
    gdf = gdf[~gdf['gene_name'].isin(genes_removed)]
    sdf = sdf[~sdf['name'].isin(genes_removed)]

    # subset the gather results to concentrate on the ones in the ground truth
    sdf_TP = sdf[sdf['name'].isin(ggene_names)]  # true positives in sourmash
    sdf_TP = sdf_TP.sort_values(by='name')  # sort by gene name
    sdf_FP = sdf[~sdf['name'].isin(ggene_names)]  # false positives in sourmash
    sdf_FN = gdf[~gdf['gene_name'].isin(sgene_names)]  # false negatives (ones in ground truth that aren't in sourmash)
    #print(f"Out of {len(sdf)} sourmash results, TP={len(sdf_TP)} are in the ground truth, FP={len(sdf_FP)} are not, "
    #      f"and there are FN={len(sdf_FN)} in the ground truth that are not in the sourmash results")
    # subset the ground truth to only the ones in the gather results
    gdf_TP = gdf[gdf['gene_name'].isin(sgene_names)]  # subset the ground truth to concentrate on the ones in the gather results
    gdf_TP = gdf_TP.sort_values(by='gene_name')
    metrics = ['TP', 'FP', 'FN', 'precision', 'recall', 'F1', 'corr_reads_mapped', 'L1_f_unique_weighted_reads_mapped', 'corr_ave_abund', 'L1_average_abund_reads_mapped_div_gene_len', 'percent_correct_predictions', 'total_number_of_predictions']
    # calculate the performance metrics
    performance = dict()
    performance['TP'] = len(sdf_TP)
    performance['FP'] = len(sdf_FP)
    performance['FN'] = len(sdf_FN)
    if performance['TP'] + performance['FP'] == 0:
        performance['precision'] = 0
    else:
        performance['precision'] = performance['TP'] / float(performance['TP'] + performance['FP'])
    if performance['TP'] + performance['FN'] == 0:
        performance['recall'] = 0
    else:
        performance['recall'] = performance['TP'] / float(performance['TP'] + performance['FN'])
    if performance['TP']:
        performance['F1'] = 2 * performance['precision'] * performance['recall'] / float(performance['precision'] + performance['recall'])
    else:
        performance['F1'] = 0
    performance['corr_reads_mapped'] = np.corrcoef(sdf_TP[sourmash_rel_abund_col], gdf_TP[ground_truth_rel_abund_col])[0][1]
    sdf_TP_vec = np.array(sdf_TP[sourmash_rel_abund_col].values)
    sdf_TP_vec = sdf_TP_vec / np.sum(sdf_TP_vec)
    gdf_TP_vec = np.array(gdf_TP[ground_truth_rel_abund_col].values)
    gdf_TP_vec = gdf_TP_vec / np.sum(gdf_TP_vec)
    performance['L1_f_unique_weighted_reads_mapped'] = np.sum(np.abs(sdf_TP_vec - gdf_TP_vec))
    reads_mapped_div_gene_len = np.array(gdf_TP['reads_mapped'] / gdf_TP['gene_length'])
    reads_mapped_div_gene_len = reads_mapped_div_gene_len / np.sum(reads_mapped_div_gene_len)
    ave_abund = sdf_TP['average_abund'] / np.sum(sdf_TP['average_abund'])
    performance['corr_ave_abund'] = np.corrcoef(ave_abund, reads_mapped_div_gene_len)[0][1]
    L1_average_abund_reads_mapped_div_gene_len = np.sum(np.abs(ave_abund - reads_mapped_div_gene_len))
    performance['L1_average_abund_reads_mapped_div_gene_len'] = L1_average_abund_reads_mapped_div_gene_len
    if len(sdf_TP) + len(sdf_FP) == 0:
        performance['percent_correct_predictions'] = 0
    else:
        performance['percent_correct_predictions'] = len(sdf_TP) / float(len(sdf_TP) + len(sdf_FP))
    performance['total_number_of_predictions'] = len(sdf_TP) + len(sdf_FP)
    # also record what filter threshold was used
    performance['filter_threshold'] = filter_threshold
    # put this in a dataframe
    performance_df = pd.DataFrame(performance, index=[0])
    return performance_df


def check_sourmash_correlation(gather_file, ground_truth_file, corr_threshold=0.9):
    """
    Since we don't know which columns of sourmash gather correlate with which columns of the ground truth, we need to
    just check them all
    :param gather_file: results of sourmash gather
    :param ground_truth_file: the output of find_genes_in_sim.py
    :param corr_threshold: only print out stats if the correlation coef is above this threshold
    :return: None
    """

    if not os.path.exists(gather_file):
        raise Exception(f"File {gather_file} does not exist")
    if not os.path.exists(ground_truth_file):
        raise Exception(f"File {ground_truth_file} does not exist")
    # s prefix is for "sourmash" while g prefix is for "ground truth"
    sdf = pd.read_csv(gather_file)
    gdf = pd.read_csv(ground_truth_file)
    # sort the ground truth by gene name, this will be the order that we stick with
    gdf = gdf.sort_values(by='gene_name')
    # grab the sourmash infered gene names
    sgene_names = list(sdf['name'])
    sgene_names = [x.split('|')[0] for x in sgene_names]
    # replace the name column with the sourmash gene names
    sdf['name'] = sgene_names
    ggene_names = list(gdf['gene_name'])
    greads_mapped = np.sum(list(gdf['reads_mapped']))
    # subset the gather results to concentrate on the ones in the ground truth
    sdf_TP = sdf[sdf['name'].isin(ggene_names)]  # true positives in sourmash
    sdf_TP = sdf_TP.sort_values(by='name')  # sort by gene name
    sdf_FP = sdf[~sdf['name'].isin(ggene_names)]  # false positives in sourmash
    sdf_FN = gdf[~gdf['gene_name'].isin(sgene_names)]  # false negatives (ones in ground truth that aren't in sourmash)
    print(f"Out of {len(sdf)} sourmash results, TP={len(sdf_TP)} are in the ground truth, FP={len(sdf_FP)} are not, "
          f"and there are FN={len(sdf_FN)} in the ground truth that are not in the sourmash results")
    # subset the ground truth to only the ones in the gather results
    gdf_TP = gdf[gdf['gene_name'].isin(sgene_names)]  # subset the ground truth to concentrate on the ones in the gather results
    gdf_TP = gdf_TP.sort_values(by='gene_name')  # sort by gene name
    # iterate over all the columns of the gather results and return the correlation coefficient with the number of reads mapped
    ground_truth_cols = ['nucleotide_overlap', 'median_coverage', 'mean_coverage', 'reads_mapped']
    sourmash_cols = ['intersect_bp', 'f_orig_query', 'f_match', 'f_unique_to_query',
                     'f_unique_weighted', 'average_abund', 'median_abund', 'std_abund',
                     'f_match_orig', 'unique_intersect_bp',
                     'gather_result_rank', 'remaining_bp']

    for gt_col in ground_truth_cols:
        for col in sourmash_cols:
            corr = np.corrcoef(sdf_TP[col], gdf_TP['reads_mapped'])[0][1]
            if corr > corr_threshold:
                print(f"gt: {gt_col}, sm:{col}: corr={corr}")

    # also look for correlation between sourmash output and #reads mapped / gene_length

    for col in sourmash_cols:
        corr = np.corrcoef(sdf_TP[col], gdf_TP['reads_mapped'] / gdf_TP['gene_length'])[0][1]
        if corr > corr_threshold:
            print(f"gt: reads mapped / gene length, sm:{col}: corr={corr}")
