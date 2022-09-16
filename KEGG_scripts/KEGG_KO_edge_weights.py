# This script will attempt to figure out branch lengths for the KO hierarchy.
# It will look at all pairwise AAI distances between descendant KO nodes of two connected nodes
# and use the median of those distances as the branch length.

import argparse
import os
import sys
from os.path import exists
import pathlib
from os import listdir
from os.path import isfile, join
import subprocess
import re
import numpy as np
from collections import Counter
import warnings
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import tempfile
import seaborn as sns
from scipy import sparse
from time import perf_counter


# get all leaf descendants of a certain node
def get_leaf_descendants(G, node):
    descendants = set()
    for n in nx.descendants(G, node):
        if G.out_degree(n) == 0:
            descendants.add(n)
    return descendants


# parse arguments
parser = argparse.ArgumentParser(description='Adds edge lengths to the KEGG hierarchy')
parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy')
parser.add_argument('-s', '--sketch_file', help='File containing all the sketches of the KOs')
parser.add_argument('-o', '--out_file', help='output file. Will be a plain text file, tsv separated by tabs.')
args = parser.parse_args()

edge_list = args.edge_list
sketch_file = args.sketch_file
out_file = args.out_file
# check that the files exist
if not exists(edge_list):
    raise FileNotFoundError(f"Could not find {edge_list}")
if not exists(sketch_file):
    raise FileNotFoundError(f"Could not find {sketch_file}")


edge_list = os.path.join("KEGG_scripts", "kegg_ko_edge_df.txt")
sketch_file = os.path.join("KEGG_scripts", "KOs_sketched_scaled_10.sig.zip")

########################################################################################################################
# Let's do the following: since I've already computed all pairwise distances, we can just make a large
# least squares problem fitting the tree distances to the pairwise distances
# Let's get the matrix describing which edges are traversed between all pairs of nodes
# This is a sparse matrix, so we'll need to use scipy.sparse

# read in the edge list
G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
G_undirected = G.to_undirected()
leaf_nodes = [node for node in G.nodes() if G.in_degree(node)!=0 and G.out_degree(node)==0]

# set the basis for the tree, which is an ordering of the edges. I'll identify each edge by its terminal node
basis = [x for x in G.nodes()]
basis_index = {node: i for i, node in enumerate(basis)}
# initialize sparse matrix on this basis
# edge_mat = sparse.lil_matrix((len(basis), len(basis)))
# Let's go with a csr_array
data = []
row_inds = []
col_inds = []

# import pairwise distances
pairwise_dist = np.load("kegg_genes_KO_1000.faa_scale_1.db.k_11_compare")
# import label names
pairwise_dist_KOs = []
with open("kegg_genes_KO_1000.faa_scale_1.db.k_11_compare.labels.txt", 'r') as f:
    for line in f.readlines():
        ko = line.strip().split('|')[-1] # KO number is the last in the pip-delim list
        ko = ko.split(':')[-1]  # remove the ko: prefix
        pairwise_dist_KOs.append(ko)
pairwise_dist_KO_index = {node: i for i, node in enumerate(pairwise_dist_KOs)}

# iterate over the pairs of nodes for which we have pairwise distances
row_ind = -1
all_paths = dict(nx.all_pairs_dijkstra_path(G_undirected))
for node_i in pairwise_dist_KOs:
    for node_j in pairwise_dist_KOs:
        row_ind += 1
        # if the nodes are the same, skip since this row of the matrix is all zeros
        if node_i != node_j:
            # get the shortest path between the two nodes
            try:
                path = all_paths[node_i][node_j]
                # get the index of each path element in the basis
                path_indices = [basis_index[node] for node in path]
                # set the corresponding entries in the sparse matrix to 1
                for path_index in path_indices:
                    data.append(1)
                    row_inds.append(row_ind)
                    col_inds.append(path_index)
            except KeyError:
                # if there is no path between the two nodes, skip
                pass
        if row_ind % 1000 == 0:
            print(f"Finished {row_ind} rows")

A = sparse.csr_matrix((data, (row_inds, col_inds)), shape=(len(pairwise_dist_KOs)**2, len(basis)))




########
# iterate over the pairs of nodes for which we have pairwise distances
tic = perf_counter()
row_ind = -1
for node_i in pairwise_dist_KOs:
    paths = nx.single_source_dijkstra_path(G_undirected, node_i)
    for node_j in pairwise_dist_KOs:
        row_ind += 1
        # if the nodes are the same, skip since this row of the matrix is all zeros
        if node_i != node_j:
            # get the shortest path between the two nodes
            try:
                path = paths[node_j]
                # get the index of each path element in the basis
                path_indices = [basis_index[node] for node in path]
                # set the corresponding entries in the sparse matrix to 1
                for path_index in path_indices:
                    data.append(1)
                    row_inds.append(row_ind)
                    col_inds.append(path_index)
            except KeyError:
                # if there is no path between the two nodes, skip
                pass
        if row_ind % 1000 == 0:
            print(f"Finished {row_ind} rows")
toc = perf_counter()
print(f"Time elapsed: {toc-tic:0.4f} seconds")  # 176.3 seconds


##################################
# Let's try to do this in parallel
