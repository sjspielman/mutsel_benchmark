# SJS
# Simulate alignments according to yeast frequencies


import sys
from pyvolve import *
from numpy import loadtxt

prefix = sys.argv[1]
freqs = loadtxt("freqs.txt")
simtree = read_tree(file = "ntaxa1024_bl0.5.tre")

partitions = []
for row in freqs:
    m = Model("mutsel", {"state_freqs": row}, scale_matrix="neutral")
    p = Partition(models = m, size = 1)
    partitions.append(p)
                        
# Simulate and save alignment in both fasta, phylip format
print "Simulating"
e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = prefix + "_simulated.fasta", ratefile = None, infofile = None)
AlignIO.convert(prefix + "_simulated.fasta", "fasta", prefix + "_simulated.phy", "phylip-relaxed")
