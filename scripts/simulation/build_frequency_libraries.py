# SJS 
# Script to extract parameters for simulation from structurally curated yeast alignments.
# Uses only alignment files with at least 150 non-gapped characters taxa.
# For assigning fitnesses to highly deleterious (i.e. not seen in alignment column) residues, we assign either a very low fitness of ~-21, or a low fitness of ~-9. These are "strongly and "weakly" deleterious, respectively. 

import os
import sys
sys.path.append("../")
import numpy as np
from Bio import AlignIO
from universal_functions import *
from pyvolve import Genetics, state_freqs
ZERO=1e-8
g = Genetics()

output_directory = "true_simulation_parameters/"
yeast_directory = "ramsey2011_yeast_alignments/" # From github repository: protein_design_and_site_variability/project_files/sequences/duncan_sequences/
taxa_threshold = 150

delfreq = {"strong":1e-9, "weak":1e-4} 
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
        
        # I am aware that this is slower.
        for delf in delfreq:

            frequencies = []
            fitnesses = []
            omegas = []    
               
            # Loop over columns, and for those with sufficient data points, derive codon frequencies and compute dN/dS
            for i in range(alnlen):
                column = str(aln[:,i])
                if (len(column) - column.count("-")) >= taxa_threshold:
                    fraw = state_freqs.ReadFrequencies("amino_acid", file = yeast_directory+file, columns=[i+1])
                    aaf = fraw.compute_frequencies()        
                    if np.count_nonzero(aaf) > 1: 
                   
                        aaf2 = np.copy(aaf)
                        aaf2[aaf2 == 0.] = delfreq[delf]
                        aaf2 /= np.sum(aaf2)
                    
                        aa_fitness = np.log(aaf2)
                        cf, cf_dict = aa_freqs_to_codon_freqs(aaf2)
                        c = dNdS_from_MutSel(cf_dict)
                        dnds = c.compute_dnds()
                               
                        frequencies.append(cf)
                        fitnesses.append(aa_fitness)
                        omegas.append(dnds)
                        
            # Save simulation information (amino acid fitnesses, codon frequencies, dN/dS)
            front = output_directory + out_prefix + "_del" + delf
            freqfile = front + "_true_codon_frequencies.txt"
            fitfile  = front + "_true_aa_fitness.txt"
            dndsfile = front + "_true_dnds.csv"
            np.savetxt(freqfile, frequencies)
            np.savetxt(fitfile, fitnesses)
            with open(dndsfile, "w") as f:
                f.write("site,dnds\n")
                for i in range(len(omegas)):
                    f.write(str(i+1)+","+str(omegas[i])+"\n")

