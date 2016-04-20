# SJS
# Simulate alignments according to yeast frequencies and deep mutational scanning data


import sys
import re
from pyvolve import *
from numpy import loadtxt
dms_mu = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}


def build_model(type, params):
    if type == "dms":
        m = Model("mutsel", {"state_freqs": params, "mu": dms_mu})
    elif type == "yeast":
        m = Model("mutsel", {"fitness": params})
    else:
        raise AssertionError("womp womp.")
    return m

name      = sys.argv[1]
paramfile = sys.argv[2]
treefile  = sys.argv[3]


params = loadtxt(paramfile)
partitions = []

if name == "HA" or name == "NP":
    type = "dms"
else:
    type = "yeast"
    
for row in params:
    m = build_model(type, row)
    p = Partition(models = m, size = 1)
    partitions.append(p)
                
simtree = read_tree(file = treefile)

bl_match = re.search("n512_(bl\d\.\d+)\.tre", treefile)
bl = bl_match.group(1)
outname = name + "_" + bl + ".phy"
e = Evolver(partitions = partitions, tree = simtree)
e(seqfile = outname, seqfmt = "phylip-relaxed", ratefile = None, infofile = None)
