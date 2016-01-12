# SJS 
# Script to extract parameters for simulation from structurally curated yeast alignments.
# Uses only alignment files with at least 150 non-gapped characters taxa.
# For assigning fitnesses to highly deleterious (i.e. not seen in alignment column) residues, we assign either a very low fitness of ~-21, or a low fitness of -6. These are "strongly and "weakly" deleterious, respectively. 
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



def save_all_info(name, frequencies, fitnesses, omegas):
    '''
        Save simulation parameters to files.
    '''
    freqfile = name + "_true_codon_frequencies.txt"
    fitfile  = name + "_true_aa_fitness.txt"
    selcfile = name + "_true_selcoeffs.csv"
    dndsfile = name + "_true_dnds.csv"

    np.savetxt(freqfile, frequencies)
    np.savetxt(fitfile, fitnesses)
    calculate_save_coeffs(fitnesses, selcfile)

    with open(dndsfile, "w") as f:
        f.write("site,dnds\n")
        for i in range(len(omegas)):
            f.write(str(i+1)+","+str(omegas[i])+"\n")



output_directory = "true_simulation_parameters/"
yeast_directory = "ramsey2011_yeast_alignments/" # From github repository: protein_design_and_site_variability/project_files/sequences/duncan_sequences/
taxa_threshold = 150

mu_dict        = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
strong_freq    = 1e-9
weak_fit       = [-6.0, -4.5] # Lowest fitness, highest fitness   
weak_threshold = -10 # Any fitness below here gets increased to an interval in weak_fit

yeastfiles = os.listdir(yeast_directory)
for file in yeastfiles:
    if file.endswith(".fasta"):
        aln = AlignIO.read(yeast_directory + file, "fasta")
        numseq = len(aln)
        alnlen = len(aln[0])
        out_prefix = file.split("_Aligned_")[0]
        
        # Skip alignment files which are too small
        if numseq < taxa_threshold:  
            continue
   
        print out_prefix

        frequencies_strong = []
        fitnesses_strong   = []
        omegas_strong      = []    
        frequencies_weak   = []
        fitnesses_weak     = []
        omegas_weak        = []  
        
        
        # Loop over columns, and for those with sufficient data points, derive codon frequencies and compute dN/dS
        for i in range(alnlen):
            column = str(aln[:,i])
            if (len(column) - column.count("-")) >= taxa_threshold:
                fraw = state_freqs.ReadFrequencies("amino_acid", file = yeast_directory+file, columns=[i+1])
                aaf = fraw.compute_frequencies()        
                if np.count_nonzero(aaf) > 1: 
               
                    #### Highly deleterious ####
                    aaf2 = np.copy(aaf)
                    aaf2[aaf2 == 0.] = strong_freq
                    aaf2 /= np.sum(aaf2)
                
                    aa_fitness = np.log(aaf2)
                    cf, cf_dict = aa_freqs_to_codon_freqs(aaf2)
                    c = dNdS_from_MutSel(cf_dict)
                    dnds = c.compute_dnds()
                           
                    frequencies_strong.append(cf)
                    fitnesses_strong.append(aa_fitness)
                    omegas_strong.append(dnds)
                    
                    #### Weakly deleterious ####
                    random_fit = np.random.uniform(low = weak_fit[0], high = weak_fit[1], size = np.sum(aafit <= weak_threshold))
                    aa_fitness[aa_fitness <= weak_threshold] = random_fit  
                    codon_freqs = codon_freqs_from_fitness_eigenvector(aa_fitness, mu_dict)
                    c = dNdS_from_MutSel(codon_freqs)
                    dnds = c.compute_dnds()
                    frequencies_weak.append(cf)
                    fitnesses_weak.append(aa_fitness)
                    omegas_weak.append(dnds)

            
            name_strong = output_directory + out_prefix + "_delstrong"
            name_weak   = output_directory + out_prefix + "_delweak"
            save_all_info(name_strong, frequencies_strong, fitnesses_strong, omegas_strong)      
            save_all_info(name_weak, frequencies_weak, fitnesses_weak, omegas_weak)      
