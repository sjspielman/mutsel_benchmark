# SJS
# Simulate alignments according to inferred empirical parameters, by swMutSel without penalty. Note that branch lengths must be appropriately re-scaled to correspond to empirical nucleotide branch lengths.
# Usage: python simulate_from_empirical.py <repository path> <dataset index>
# NOTE: universal_functions.py must be in working directory!

import os
import sys
import numpy as np
from pyvolve import *
from dendropy import Tree
from universal_functions import *

datasets = ["PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715"] 

# Input arguments
repodir = sys.argv[1]
assert(os.path.exists(repodir)), "\nIncorrect repository path specified."
if not repodir.endswith("/"):
    repodir += "/"
index = int(sys.argv[2]) - 1
assert(index >= 0 and index < len(datasets)), "\nIncorrect dataset index provided."


# Directories
outdir  = repodir + "data/simulated_from_empirical/"
resdir  = repodir + "results/raw_results/empirical/swmutsel/"
treedir = repodir + "data/empirical/"


# Dataset for simulation
dataset = datasets[index]

# Inferred information about dataset
fitness  = np.loadtxt(resdir + dataset + "_nopenal_fitness.txt")
mu_dict, scaling  = parse_swMutSel_mutation(resdir + dataset + "_nopenal_MLE.txt", return_scaling = True)   

# Update tree with appropriately scaled branch lengths, remove unrooted silliness tag, and save
treefile_raw = treedir + dataset + ".tre"
treefile  = outdir + dataset + "_scaled.tre"
rawtree = Tree.get_from_path(treefile_raw, "newick")
rawtree.scale_edges(scaling)
treestring = str(rawtree).replace('[&U] ','') + ";"
with open(treefile, "w") as f:
    f.write(treestring)			

# Simulate with pyvolve
simtree = read_tree(tree = treestring)
partitions = []
for row in fitness:
    m = Model("mutsel", {"fitness": row, "mu": mu_dict})
    p = Partition(models = m, size = 1)
    partitions.append(p)
                        
outfile1 = outdir + dataset + "_simulated.fasta"
outfile2 = outdir + dataset + "_simulated.phy"

e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = outfile1, ratefile = None, infofile = None)
AlignIO.convert(outfile1, "fasta", outfile2, "phylip-relaxed")





