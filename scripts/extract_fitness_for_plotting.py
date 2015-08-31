# SJS
# Script to extract state frequencies from fitness values across sites for a representative simulated dataset to make a comparative barplot.

from numpy import *
from parsing_functions import *
import sys

usage = "Usage: python extract_fitness_for_plotting.py <dataset> , where dataset is the name of a simulated dataset."
assert(len(sys.argv) == 2), usage



dataset = sys.argv[1]
dataset = dataset.replace("_simulated", "")

outfile = "postprocessing/fitness_comparison_" + dataset + "_simulated.csv" 
# True frequencies
true_freqs_file = "simulation/flib/" + dataset + "_codon_freq_lib.txt"

# swMutSel raw data
nopenal_fitfile = "../results/simulation/" + dataset + "_simulated_nopenal_fitness.txt"

# pbMutSel raw data
#pb_fitfile      = "../results/simulation/" + dataset + "_phylobayes.aap"

# Parse swMutSel and determine the relative value to set to 0
print "Parsing swMutSel"

nopenal_fitness = np.loadtxt(nopenal_fitfile)

# Unfortunately, have to do this in a loop.
nopenal_fitness_norm = []
#nopenal_fitness_norm = (nopenal_fitness.T / nopenal_fitness.max(axis=1)).T
zero_fitness_indices = []
x = 0

for row in nopenal_fitness:
    
    # Zero fitness entry
    ind = np.where(row == 0.)[0]
    zero_fitness_indices.append( ind )
    
    norm_row = None
    row += np.abs(np.min(row))
    # Normalization
    if np.max(row) == 0.:
        row2 = row + 5
        norm_row = row2/np.max(row2)
    else:
        #print row
        norm_row = row/np.max(row)

    assert(norm_row is not None), "couldn't normalize row" 
    nopenal_fitness_norm.append(norm_row)


# Turn true codon frequencies into amino acid relative fitness values with same "zero fitness" as swmutsel
print "Parsing true"
true_fitness_norm = []
true_cfreqs = np.loadtxt(true_freqs_file)
i = 0
for row in true_cfreqs:
    row[row <= ZERO] = 1e-15 # replace zeros with tiny values, since log on next line
    cfit = np.log(row)
    
    aafit = []
    for syn in g.genetic_code:
        synindex = g.codons.index(syn[0])
        aafit.append( cfit[synindex] )
    aafit += np.abs(np.min(aafit))
    factor = 0. - aafit[ zero_fitness_indices[i] ]
    aafit -= factor
    aafit_norm = aafit/np.max(aafit)
    true_fitness_norm.append(aafit_norm)

    

print "Saving"
with open(outfile, "w") as f:
    f.write("site,method,fitness,aminoacid\n")
    for i in range(len(true_fitness_norm)):
        swmutsel_row = nopenal_fitness_norm[i]
        true_row     = true_fitness_norm[i]
        for j in range(len(g.amino_acids)):
        
            line1 = str(i) + ",nopenal," + str(swmutsel_row[j]) + "," + g.amino_acids[j] + "\n"
            line2 = str(i) + ",true," + str(true_row[j]) + "," + g.amino_acids[j] + "\n"
            f.write(line1 + line2)
            

           
# # R code for barplot            
# for (i in 1:nrow(dat)){
# dat %>% filter(site == i) %>% ggplot(aes(x = aminoacid, y = fitness, group = method, fill = method)) + geom_bar(position="dodge", stat="identity", width = 0.5)  + ggtitle(paste0("site", i)) -> p
# print(p)
# readline("enter to continue")
# }    
# 
      

















