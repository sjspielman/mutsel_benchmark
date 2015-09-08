# SJS
# Here, we predict site-specific dN/dS values from the frequencies in the alignment columns.

from compute_dnds_from_mutsel import *
from universal_functions import *
from pyvolve import Genetics
from Bio import AlignIO
import numpy as np
import os
ZERO=1e-10
g = Genetics()




def compute_dnds(d, mu_dict):
    c = dNdS_from_MutSel(d, mu_dict)
    dnds = None
    # Assertion error will be raised in dnds calculation if no evolution has occurred. These are uninformative sites.
    try:
        dnds = c.compute_dnds()
        if np.isnan(dnds) or np.isinf(dnds):
            dnds = "NA"
    except AssertionError:
        dnds = "NA"
    
    assert(dnds != None), "dnds not computed."
    return dnds


def derive_site_dnds(site, records, mu_dict):
    
    freq_dict = dict(zip(g.codons, np.zeros(61)))
    freq_dict_from_aa = {}  # dict(zip(codons, np.repeat(0., 61)))
    
    # Collect codon frequencies, ignoring any gaps or noncanonical codons that may exist. Throw GIANT error upon encountering a stop codon
    column = []
    print site
    column_raw = list(records[:, site:site+3])
    for entry in column_raw:
        codon_current = str(entry.seq)
        assert(codon_current not in g.stop_codons), "WOAH, STOP CODON!"
        if codon_current in g.codons:
            column.append(codon_current)
    length = len(column)
    for codon in column:
        freq_dict[codon] += 1./length
    
    # Determine corresponding amino acid frequencies and re-convert to codon. This removes signal of synonymous variation from alignment for better comparison with swMutSel, pbMutSel
    for aa_family in g.genetic_code:
        total_aa = 0.
        for entry in aa_family:
            total_aa += freq_dict[entry]
        total_aa /= float(len(aa_family))
        for entry in aa_family:
            freq_dict_from_aa[entry] = total_aa
         
    # Calculate the dN/dS  
    dnds = compute_dnds(freq_dict_from_aa, mu_dict)       

    return dnds


emp_datasets = ["PF00126", "PF04055", "PF00593","PF07715", "PF02518", "PF01336", "PF01926","PF01266",  "pb2" ,"amine"]
sim_datasets = ["1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
datasets = {"simulation":sim_datasets}

for datatype in datasets:
    
    seqdir = "../data/" + datatype + "/"
    resdir = "../results/" + datatype + "/"
    
    for data in datasets[datatype]:
        print data
        
        outfile = resdir + data + "_dnds_from_alignment.txt"
        alnfile = seqdir + data + ".fasta"
        hyphile = resdir + data + "_swmutsel_hyout.txt"

        if datatype == "empirical":
            mu_dict = extract_optimized_params(hyphile, return_tree = False)
        else:
            mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
    
        records = AlignIO.read(alnfile, "fasta")
        omegas = []
    
        for i in range(len(records[0])/3):
            site = i*3
            omegas.append( derive_site_dnds(site, records, mu_dict) )
    
        
        with open(outfile, "w") as outf:
            outf.write("site,empdnds\n")
            for i in range(len(omegas)):
                outf.write(str(i) + "," + str(omegas[i]) + "\n")
    
                
                        

