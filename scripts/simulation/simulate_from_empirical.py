# SJS
# Simulate alignments according to inferred empirical parameters, by swMutSel without penalty.
# Usage: python simulate_from_inferred.py <repository path> <dataset index>
# NOTE: universal_functions.py must be in working directory!

import os
import sys
import numpy as np
from pyvolve import *
from universal_functions import *

datasets = ["PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715", "pb2"] 

# Input arguments
repodir = sys.argv[1]
assert(os.path.exists(repodir)), "\nIncorrect repository path specified."
if not repodir.endswith("/"):
    repodir += "/"
index = int(sys.argv[2]) - 1
assert(index >= 0 and index < len(datasets)), "\nIncorrect dataset index provided."


# Dataset for simulation
dataset = datasets[index]

# Inferred information about dataset
treedir  = repodir + "data/empirical/"
resdir   = repodir + "results/raw_results/empirical/swmutsel/"
treefile = treedir + dataset + ".tre"
fitness  = np.loadtxt(resdir + dataset + "_nopenal_fitness.txt")
mu_dict  = parse_swMutSel_mutation(resdir + dataset + "_nopenal_MLE.txt")   


# Simulate with pyvolve
simtree = read_tree(file = treefile)
partitions = []
for row in fitness:
    m = Model("mutsel", {"fitness": row, "mu": mu_dict})
    p = Partition(models = m, size = 1)
    partitions.append(p)
                        
outdir = repodir + "data/simulated_from_empirical/"
outfile1 = outdir + dataset + "_simulated.fasta"
outfile2 = outdir + dataset + "_simulated.phy"

e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = outfile1, ratefile = None, infofile = None)
AlignIO.convert(outfile1, "fasta", outfile2, "phylip-relaxed")





