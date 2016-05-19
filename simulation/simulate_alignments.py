"""
    SJS
    Simulate alignments according to yeast frequencies or deep mutational scanning data.
    Usage: python simulate_alignments.py <name> <freqfile> <treefile>. Name is the dataset name, freqfile is a file of codon frequencies, treefile is the file with the tree for simulation.
"""

import re
import sys
from pyvolve import *
from numpy import loadtxt
mudict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}


name     = sys.argv[1]
freqfile = sys.argv[2]
treefile = sys.argv[3]


frequencies = loadtxt(freqfile)
partitions = []

    
for row in frequencies:
    m = Model("mutsel", {"state_freqs": row, "mu": mudict})
    p = Partition(models = m, size = 1)
    partitions.append(p)
                
simtree = read_tree(file = treefile)

bl_match = re.search("n512_(bl\d\.\d+)\.tre", treefile)
bl = bl_match.group(1)
outname = name + "_" + bl + ".phy"
e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = outname, seqfmt = "phylip-relaxed", ratefile = None, infofile = None)
