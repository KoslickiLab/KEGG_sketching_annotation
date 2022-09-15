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

# read in the edge list
G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
G_undirected = G.to_undirected()

# get all leaf descendants of a certain node
def get_leaf_descendants(G, node):
    descendants = set()
    for n in nx.descendants(G, node):
        if G.out_degree(n) == 0:
            descendants.add(n)
    return descendants

# Just prototype code for now
# open a temporary file to make the picklist
#KOs_1 = get_leaf_descendants(G, "00680 Methane metabolism")
KOs_1 = get_leaf_descendants(G, "04147 Exosome")
KOs_2 = get_leaf_descendants(G, "09101 Carbohydrate metabolism")
# Picklists don't play nicely in windows due to the ':' delimiter and in file paths too
#pick_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".csv")
pick_file = "KO1_KO2_picklist.csv"
with open(pick_file, 'w') as picklist_file:
    picklist_file.write("name\n")
    for ko in KOs_1:
        picklist_file.write(f"ko:{ko}\n")
    for ko in KOs_2:
        picklist_file.write(f"ko:{ko}\n")

# run sourmash compare
sourmash_compare_cmd = f"sourmash compare -k 5 --protein --ani --picklist {pick_file}:name:name {sketch_file} -o {pick_file}.compare"
subprocess.run(sourmash_compare_cmd, shell=True)

# import the results
dist_mat = np.load(f"{pick_file}.compare")

# get the median distance between the two sets of KOs, ignoring any zeros
dist_mat_pairs = dist_mat[0:len(KOs_1), len(KOs_1):]
dist_mat_pairs_nonzero = dist_mat[dist_mat > 0]
median_dist = np.median(dist_mat_pairs_nonzero)
print(f"Median distance between KOs: {median_dist}")
sns.distplot(dist_mat_pairs_nonzero)
plt.show()

########################################################################################################################
# Alternate approach: union up the KO sketches, then compare them
#KOs_1 = get_leaf_descendants(G, "00680 Methane metabolism")
KOs_1 = get_leaf_descendants(G, "04147 Exosome")
KOs_2 = get_leaf_descendants(G, "09101 Carbohydrate metabolism")
#pick_file_1 = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".csv")
#pick_file_2 = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".csv")
pick_file_1 = "KO1_picklist.csv"
pick_file_2 = "KO2_picklist.csv"
with open(pick_file_1, 'w') as picklist_file:
    picklist_file.write("name\n")
    for ko in KOs_1:
        picklist_file.write(f"ko:{ko}\n")
with open(pick_file_2, 'w') as picklist_file:
    picklist_file.write("name\n")
    for ko in KOs_2:
        picklist_file.write(f"ko:{ko}\n")
# merge the picklists
cmd = f"sourmash sig merge -o {pick_file_1}.merged.sig --protein -k 5 --picklist {pick_file_1}:name:name {sketch_file}"
subprocess.run(cmd, shell=True)
cmd = f"sourmash sig merge -o {pick_file_2}.merged.sig --protein -k 5 --picklist {pick_file_2}:name:name {sketch_file}"
subprocess.run(cmd, shell=True)
# run sourmash compare
csv_file = f"{pick_file_1}_{pick_file_2}.compare.csv"
cmd = f"sourmash compare -k 5 --protein --ani {pick_file_1}.merged.sig {pick_file_2}.merged.sig --csv {csv_file}"
res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
# import the matrix
dist_mat = np.loadtxt(csv_file, delimiter=',', skiprows=1)
dist = dist_mat[0, 1]
print(f"Distance between union of KOs: {dist}")




########################################################################################################################
# Let's do the following: since I've already computed all pairwise distances, we can just make a large
# least squares problem fitting the tree distances to the pairwise distances
# Let's get the matrix describing which edges are traversed between all pairs of nodes
# This is a sparse matrix, so we'll need to use scipy.sparse

edge_list = os.path.join("KEGG_scripts", "kegg_ko_edge_df.txt")
sketch_file = os.path.join("KEGG_scripts", "KOs_sketched_scaled_10.sig.zip")

# read in the edge list
G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
G_undirected = G.to_undirected()
leaf_nodes = [node for node in G.nodes() if G.in_degree(node)!=0 and G.out_degree(node)==0]
KOs_list_small = leaf_nodes[:100]
# set the basis for the tree, which is an ordering of the edges. I'll identify each edge by its terminal node
basis = [x for x in G.nodes()]
basis_index = {node: i for i, node in enumerate(basis)}
# initialize sparse matrix on this basis
edge_mat = sparse.lil_matrix((len(basis), len(basis)))

# iterate over the pairs of nodes
for node_i in KOs_list_small:
    for node_j in KOs_list_small:
        # get the shortest path between the two nodes
        path = nx.shortest_path(G_undirected, node_i, node_j)
        # get the index of each path element in the basis
        path_indices = [basis_index[node] for node in path]
        # set the corresponding entries in the sparse matrix to 1
        edge_mat[path_indices, path_indices] = 1

