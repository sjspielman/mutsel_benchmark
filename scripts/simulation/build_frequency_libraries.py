# Grabs amino acid frequencies from the Ramsey paper. Uses only alignment files with at least 100 non-gapped characters taxa.

import os
import sys
sys.path.append("../") #for compute_dnds_from_mutsel module
import subprocess
from Bio import AlignIO
from numpy import savetxt, count_nonzero
from compute_dnds_from_mutsel import *
from pyvolve import Genetics, state_freqs
ZERO=1e-8
g = Genetics()

yeast_directory = "ramsey2011_yeast_alignments/" # From github repository: protein_design_and_site_variability/project_files/sequences/duncan_sequences/
taxa_threshold = 100

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
        frequencies = []
        omegas = []   
        
        # Loop over columns, and for those with sufficient data points, derive codon frequencies and compute dN/dS
        for i in range(alnlen):
            column = str(aln[:,i])
            if (len(column) - column.count("-")) >= taxa_threshold:
                fraw = state_freqs.ReadFrequencies("amino_acid", file = yeast_directory+file, columns=[i+1])
                aaf = fraw.compute_frequencies()
                cf = fraw.compute_frequencies(type = "codon")
                
                # Only save the frequencies if contains multiple amino acids, the different codons are accessible through *single* instantaneous changes, and dN/dS is neither 0 nor NA
                if count_nonzero(aaf) > 1:
                    cf_dict = dict(zip(g.codons, cf))
                    save = True
                    try:    
                        c = dNdS_from_MutSel(cf_dict)
                        dnds = c.compute_dnds()
                    except:
                        save = False

                    # Retain this information
                    if save and dnds > ZERO:
                        #saving
                        frequencies.append(cf)
                        omegas.append(dnds)
            
        # Save simulation information (true omegas and frequencies)
        freqfile = out_prefix + "_codon_freq_lib.txt"
        dndsfile = out_prefix + "_true_dnds.txt"
        savetxt(freqfile, frequencies)
        with open(dndsfile, "w") as f:
            f.write("site,dnds\n")
            for i in range(len(omegas)):
                f.write(str(i+1)+","+str(omegas[i])+"\n")
            
