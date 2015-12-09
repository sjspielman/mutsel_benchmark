# SJS
# Script to show how a given set of equilibrium frequencies may have different fitness distributions, and hence different distribitions of selection coefficients, depending on how deleterious the deleterious changes are. Studying these distributions is counterproductive for understanding site-specific evolutionary constraint.

import numpy as np
import sys
sys.path.append("../")
from universal_functions import *

def calc_entropy(f):
    return -1. * np.sum ( f[f > 1e-8] * np.log(f[f > 1e-8]) )    


data = "1B4T_A"
flib = "simulation/flib/" + data + "_codon_freq_lib.txt"
original_freqs = np.loadtxt(flib)

kappa=1.;a=1;c=1.;g=1.;t=1.;
mu_dict = {'AG':kappa*g, 'TC':kappa*c, 'GA':kappa*a, 'CT':kappa*t, 'AC':c, 'TG':g, 'CA':a, 'GT':t, 'AT':t, 'TA':a, 'GC':c, 'CG':g}  

outfile = data + "_vary_deleterious.txt"

site = 1
with open(outfile, "w") as f:
    f.write("site,scalar,dnds\n")
    
for origf in original_freqs:
    print site
    fitness = np.log(origf)
    min_fitness = np.min(fitness[np.isfinite(fitness)])
    masker = np.isinf(fitness)
    
    c = dNdS_from_MutSel(origf)
    true_dnds = c.compute_dnds()
    with open(outfile, "a") as f:
        f.write(str(site) + ",true," + str(true_dnds) + "\n")
        
    for scalar in np.arange(1.5, 10, 0.5):
        temp = np.copy(fitness)
        temp[masker] = scalar * min_fitness
        dnds = dnds_from_params(temp, mu_dict)
        with open(outfile, "a") as f:
            f.write(str(site) + "," + str(scalar) + "," + str(dnds) + "\n")
    
    site += 1














