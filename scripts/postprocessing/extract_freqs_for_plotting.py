# SJS
# Script to extract state frequencies from fitness values across sites for a representative simulated dataset to make a comparative barplot.

import sys
sys.path.append("../") # universal_functions.py is one directory up
import numpy as np
from universal_functions import *


dataset = "1IBS_A" # Same representative simulation set used for scatterplots
outfile = "aa_frequency_comparison_" + dataset + "_simulated.csv" 

# True frequencies, dN/dS
true_freqs_file = "../simulation/flib/" + dataset + "_codon_freq_lib.txt"
true_dnds_file = "../../results/simulation/" + dataset + "_true_dnds.txt"

# swMutSel raw data
nopenal_fitfile = "../../results/simulation/" + dataset + "_simulated_nopenal_fitness.txt"
nopenal_mlefile = "../../results/simulation/" + dataset + "_simulated_nopenal_MLE.txt"

# pbMutSel raw data
pb_fitfile      = "../../results/simulation/" + dataset + "_simulated_phylobayes.aap"
pb_tracefile    = "../../results/simulation/" + dataset + "_simulated_phylobayes.trace"


# Obtain true dN/dS
with open(true_dnds_file) as f:
    truednds_raw = f.read().strip().split("\n")
truednds = []
for raw in truednds_raw[1:]:
    truednds.append(raw.split(",")[1])

# Parse swMutSel
print "Parsing swMutSel"
swmutsel_aafreqs = []
sitewise_fitness = np.loadtxt(nopenal_fitfile)
mu_dict = parse_swMutSel_mutation(nopenal_mlefile)
for site_fitness in sitewise_fitness:
      eqfreqs = codon_freqs_from_fitness_eigenvector(site_fitness, mu_dict)
      aafreqs = codon_freqs_to_aa_freqs(eqfreqs)
      swmutsel_aafreqs.append(aafreqs)

# Parse pbMutSel
print "Parsing pbMutSel"
pbmutsel_aafreqs = []
sitewise_fitness = np.log( np.loadtxt(pb_fitfile) )
mu_dict = parse_pbMutSel_mutation(pb_tracefile)    
for site_fitness in sitewise_fitness:
      eqfreqs = codon_freqs_from_fitness_eigenvector(site_fitness, mu_dict)
      aafreqs = codon_freqs_to_aa_freqs(eqfreqs)
      pbmutsel_aafreqs.append(aafreqs)


# Turn true codon frequencies into amino acid frequencies
print "Parsing true"
true_aafreqs = []
true_cfreqs = np.loadtxt(true_freqs_file)
for row in true_cfreqs:
    true_aafreqs.append( codon_freqs_to_aa_freqs(row) )

# Check size
assert(len(true_aafreqs) == len(swmutsel_aafreqs)), "Different number of positions for true, swmutsel..."
assert(len(true_aafreqs) == len(pbmutsel_aafreqs)), "Different number of positions for true, pnmutsel..."


print "Saving"
with open(outfile, "w") as f:
    f.write("site,truednds,method,freq,aminoacid\n")
    for i in range(len(true_aafreqs)):
        true_row     = true_aafreqs[i]
        swmutsel_row = swmutsel_aafreqs[i]
        pbmutsel_row = pbmutsel_aafreqs[i]
        for j in range(len(g.amino_acids)):
        
            line1 = str(i+1) + "," + truednds[i] + ",true," + str(true_row[j]) + "," + g.amino_acids[j] + "\n"
            line2 = str(i+1) + "," + truednds[i] + ",swmutsel," + str(swmutsel_row[j]) + "," + g.amino_acids[j] + "\n"
            line3 = str(i+1) + "," + truednds[i] + ",pbmutsel," + str(pbmutsel_row[j]) + "," + g.amino_acids[j] + "\n"

            f.write(line1 + line2 + line3)
            
            
            
    

      

















