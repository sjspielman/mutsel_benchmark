# SJS
# Simulate alignments according to yeast frequencies


import sys
from pyvolve import *
from numpy import loadtxt


name     = sys.argv[1]
fitfile  = sys.argv[2]
treefile = sys.argv[3]

aa_fitnesses = loadtxt(fitfile)


partitions = []
for row in aa_fitnesses:
    m = Model("mutsel", {"fitness": row})
    p = Partition(models = m, size = 1)
    partitions.append(p)
                
simtree = read_tree(file = treefile)
outname = name + "_simulated.phy"
e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = outname + ".phy", seqfmt = "phylip-relaxed", ratefile = None, infofile = None)
