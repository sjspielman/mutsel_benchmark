# SJS 
# Script to extract parameters for simulation from structurally curated yeast alignments.
# Uses only alignment files with at least 150 non-gapped characters taxa.
# For assigning fitnesses to highly deleterious (i.e. not seen in alignment column) residues, we assign either a very low fitness of ~-21, or a low fitness in the [-6,-4.5] interval. These are "strongly and "weakly" deleterious, respectively. 
# Saves site-specific amino acid fitnesses, site-specific equilibrium codon frequencies, site-specific true dN/dS, and selection coefficient distributions

import os
import sys
sys.path.append("../")
import numpy as np
from Bio import AlignIO
from universal_functions import *
from pyvolve import Genetics, state_freqs
ZERO=1e-8
g = Genetics()

def calculate_save_coeffs(fitness, outfile):
    '''
        Compute and save distribution of selection coefficients.
    '''
    raw = []
    binned = []
    for site in fitness:
        for i in range(len(site)):
            f_i = site[i]
            for j in range(len(site)):
                if i == j:
                    continue
                else:
                    s = f_i - site[j]
                    raw.append(s)
                    if s <= -10:
                        s = -10.
                    if s >= 10.:
                        s = 10.
                    binned.append(s)
    with open(outfile, "w") as outf:
        outf.write("realcoeff,binnedcoeff\n")
        outf.write( "\n".join([str(raw[x])+","+str(binned[x]) for x in range(len(raw))]) )

output_directory = "true_simulation_parameters/"
names = ["1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A"]
mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}

low_weak = -6
high_weak = -4.5

for name in names:
    print name
    omegas = []
    frequencies = []
    fitnesses = []
    fitfile = output_directory + name + "_delweak_true_aa_fitness.txt"
    freqfile = output_directory + name + "_delweak_true_codon_frequencies.txt"
    dndsfile = output_directory + name + "_delweak_true_dnds.csv"
    selcfile = output_directory + name + "_delweak_true_selcoeffs.csv"
    
    origfit = np.loadtxt( output_directory + name + "_delstrong_true_aa_fitness.txt" )
    for row in origfit:
    
        aafit = np.copy(row)
        random_fit = np.random.uniform(low = -6, high = -4.5, size = np.sum(aafit <= -10))
        aafit[aafit <= -10] = random_fit  
        
                    
        codon_freqs = codon_freqs_from_fitness_eigenvector(aafit, mu_dict)
        c = dNdS_from_MutSel(codon_freqs)
        dnds = c.compute_dnds()

        omegas.append(dnds)
        frequencies.append(codon_freqs)
        fitnesses.append(aafit)
    
    np.savetxt(freqfile, frequencies)
    np.savetxt(fitfile, fitnesses)
    with open(dndsfile, "w") as f:
        f.write("site,dnds\n")
        for i in range(len(omegas)):
            f.write(str(i+1)+","+str(omegas[i])+"\n")
    calculate_save_coeffs(fitnesses, selcfile)

