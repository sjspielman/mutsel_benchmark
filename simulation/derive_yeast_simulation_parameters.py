# SJS 
# Script to extract parameters for simulation from structurally curated yeast alignments.
# Uses only alignment files with at least 150 non-gapped characters taxa.
# For assigning fitnesses to highly deleterious (i.e. not seen in alignment column) residues, we assign either a very low fitness of ~-21, or a low fitness of -6. These are "strongly and "weakly" deleterious, respectively. 
# Saves site-specific amino acid fitnesses, site-specific equilibrium codon frequencies, site-specific true dN/dS, and selection coefficient distributions

import os
import sys
import numpy as np
from Bio import AlignIO
from pyvolve import Genetics, state_freqs

sys.path.append("../")
from universal_functions import *


output_directory = "true_simulation_parameters/"
yeast_directory = "ramsey2011_yeast_alignments/" # From github repository: protein_design_and_site_variability/project_files/sequences/duncan_sequences/
taxa_threshold = 150

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
                    aaf2[aaf2 == 0.] = STRONG_FREQ
                    aaf2 /= np.sum(aaf2)
                
                    aa_fitness = np.log(aaf2)
                    cf, cf_dict = aa_freqs_to_codon_freqs(aaf2)
                    c = dNdS_from_MutSel(cf_dict)
                    dnds = c.compute_dnds()
                           
                    frequencies_strong.append(cf)
                    fitnesses_strong.append(aa_fitness)
                    omegas_strong.append(dnds)
                    
                    #### Weakly deleterious ####
                    random_fit = np.random.uniform(low = WEAK_FIT[0], high = WEAK_FIT[1], size = np.sum(aa_fitness <= WEAK_THRESHOLD))
                    aa_fitness[aa_fitness <= WEAK_THRESHOLD] = random_fit  
                    codon_freqs = codon_freqs_from_fitness_eigenvector(aa_fitness, mu_dict)
                    c = dNdS_from_MutSel(codon_freqs)
                    dnds = c.compute_dnds()
                    frequencies_weak.append(cf)
                    fitnesses_weak.append(aa_fitness)
                    omegas_weak.append(dnds)

            
            name_strong = output_directory + out_prefix + "_delstrong"
            name_weak   = output_directory + out_prefix + "_delweak"
            save_simulation_info(name_strong, frequencies_strong, fitnesses_strong, omegas_strong)      
            save_simulation_info(name_weak, frequencies_weak, fitnesses_weak, omegas_weak)      
