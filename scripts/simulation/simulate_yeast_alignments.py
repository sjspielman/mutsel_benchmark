# SJS
# Simulate alignments according to yeast frequencies


import sys
from pyvolve import *
from numpy import loadtxt

prefix = sys.argv[1]
freqs = loadtxt("freqs.txt")
simtree = read_tree(file = "ntaxa256_bl0.5.tre")

partitions = []
for row in freqs:
    m = Model("mutsel", {"state_freqs": row})
    p = Partition(models = m, size = 1)
    partitions.append(p)
                        
# Simulate and save alignment in both fasta, phylip format
print "Simulating"
e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = prefix + "_simulated_n256.fasta", ratefile = None, infofile = None)
AlignIO.convert(prefix + "_simulated_n256.fasta", "fasta", prefix + "_simulated_n256.phy", "phylip-relaxed")
