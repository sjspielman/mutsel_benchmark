"""
    SJS
    Script to extract parameters for simulation from deep mutational scanning data. Some notes..
        Extract codon frequencies, dnds from deep mutational scanning AA preferences combined with experimental mutation rates. 
        Note that these parameters, as they are experimental, are simply taken as they are, and NOT subjected to a strongly/weakly deleterious regime.
"""

import sys
sys.path.append("../")
from universal_functions import *
import numpy as np




mudict  = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
truedir = "true_simulation_parameters/"

for source in ["HA", "NP", "LAC", "Gal4"]:
    print source
    infile  = truedir + source + "_prefs.txt"   
    outname = truedir + source   
    rawpref = np.loadtxt(infile)
    nsites  = len(rawpref)

    final_codon_freqs = np.zeros([nsites, 61])
    final_fitness = np.zeros([nsites, 20])
    final_dnds = np.zeros(nsites)
    final_entropy = np.zeros(nsites)
    for i in range(nsites):
        print "  ",i
        aa_freqs      = rawpref[i] / np.sum(rawpref[i]) # Renormalize. They are at a bad tolerance.
        aa_fitness    = np.log(aa_freqs)
        codon_fitness = aa_fitness_to_codon_fitness(aa_fitness)
        cf = codon_freqs_from_fitness_boltzmann(codon_fitness)
        c = dNdS_from_MutSel(dict(zip(g.codons,cf)), mudict)
        dnds = c.compute_dnds()

        final_codon_freqs[i] = cf
        final_fitness[i]     = aa_fitness
        final_dnds[i]        = dnds
        final_entropy[i]     = calculate_entropy( aa_freqs ) 


    save_simulation_info(outname, final_codon_freqs, final_fitness, final_dnds, final_entropy)
