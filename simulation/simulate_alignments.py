"""
    SJS
    Simulate alignments according to yeast frequencies or deep mutational scanning data.
    Usage: python simulate_alignments.py <name> <freqfile> <treefile>. Name is the dataset name, freqfile is a file of codon frequencies, treefile is the file with the tree for simulation.
"""

import re
import sys
from pyvolve import *
from numpy import loadtxt
h3n2_mudict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
equal_mudict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}


name     = sys.argv[1]
freqfile = sys.argv[2]
treefile = sys.argv[3]


frequencies = loadtxt(freqfile)
partitions = []


if name in ["HA", "NP"]:
    mudict = h3n2_mudict
else:
    mudict = equal_mudict
    
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
