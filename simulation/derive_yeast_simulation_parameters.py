"""
    SJS 
    Script to extract parameters for simulation from Ramsey et al. (2011) yeast alignments. Some notes..
        Uses only alignment files with at least 150 non-gapped characters taxa.
        For assigning fitnesses to highly deleterious (i.e. not seen in alignment column) residues, we assign either a very low fitness of ~-21, or a low fitness of -6. These are "strongly and "weakly" deleterious, respectively. 
        Saves site-specific amino acid fitnesses, site-specific equilibrium codon frequencies, site-specific true dN/dS, and selection coefficient distributions
"""

import os
import sys
import numpy as np
from Bio import AlignIO
from pyvolve import state_freqs
sys.path.append("../")
from universal_functions import *

### Simulation constants and functions ###
mu_dict        = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
STRONG_FREQ    = 1e-9

def apply_weakdel(fitness):
    '''
        Convert a fitness distribution to weakly deleterious regime.
    ''' 
    random_fit = np.random.uniform(low = WEAK_FIT[0], high = WEAK_FIT[1], size = np.sum(fitness <= WEAK_THRESHOLD))
    fitness[fitness <= WEAK_THRESHOLD] = random_fit  
    return fitness
            
            
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
        outname = output_directory + out_prefix

        frequencies = []
        fitnesses   = []
        omegas      = []   
        entropies   = [] 
        
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
                    c = dNdS_from_MutSel(dict(zip(g.codons, cf)), mu_dict)
                    dnds = c.compute_dnds()
                           
                    frequencies.append(cf)
                    fitnesses.append(aa_fitness)
                    omegas.append(dnds)
                    entropies.append(calculate_entropy( codon_freqs_to_aa_freqs(cf)) )
             
        save_simulation_info(outname, frequencies, fitnesses, omegas, entropies)      
