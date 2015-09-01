# SJS
# Script to extract state frequencies from fitness values across sites for a representative simulated dataset to make a comparative barplot.

from numpy import *
from parsing_functions import *
import sys


dataset = "1HUR_A"
outfile = "postprocessing/aa_frequency_comparison_" + dataset + "_simulated.csv" 
# True frequencies
true_freqs_file = "simulation/flib/" + dataset + "_codon_freq_lib.txt"

# swMutSel raw data
hyphy_file      = "../results/simulation/" + dataset + "_simulated_swmutsel_hyout.txt"
nopenal_fitfile = "../results/simulation/" + dataset + "_simulated_nopenal_fitness.txt"

# pbMutSel raw data
#pb_fitfile      = "../results/simulation/" + dataset + "_phylobayes.aap"
#pb_tracefile    = "../results/simulation/" + dataset + "_phylobayes.trace"

# Parse swMutSel
print "Parsing swMutSel"
swmutsel_aafreqs = []
sitewise_fitness = np.loadtxt(nopenal_fitfile)
pi, kappa, mu_dict = extract_optimized_params(hyphy_file, return_tree = False)       
for site_fitness in sitewise_fitness:
      eqfreqs = codon_freqs_from_fitness(site_fitness, mu_dict)
      aafreqs = codon_to_aa_freqs(eqfreqs)
      swmutsel_aafreqs.append(aafreqs)

# Turn true codon frequencies into amino acid frequencies
print "Parsing true"
true_aafreqs = []
true_cfreqs = np.loadtxt(true_freqs_file)
for row in true_cfreqs:
    true_aafreqs.append( codon_to_aa_freqs(row) )

# Check size
assert(len(true_aafreqs) == len(swmutsel_aafreqs)), "Different number of positions for true, swmutsel..."

print "Saving"
with open(outfile, "w") as f:
    f.write("site,method,freq,aminoacid\n")
    for i in range(len(true_aafreqs)):
        swmutsel_row = swmutsel_aafreqs[i]
        true_row     = true_aafreqs[i]
        for j in range(len(g.amino_acids)):
        
            line1 = str(i) + ",nopenal," + str(swmutsel_row[j]) + "," + g.amino_acids[j] + "\n"
            line2 = str(i) + ",true," + str(true_row[j]) + "," + g.amino_acids[j] + "\n"
            f.write(line1 + line2)
            
            
            
    

      

















