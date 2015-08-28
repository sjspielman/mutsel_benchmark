# SJS
# Here, we predict site-specific dN/dS values from the frequencies in the alignment columns.
# We need to use the mutation rates as inferred by HyPhy for this, as well.

from compute_dnds_from_mutsel import *
from pyvolve import Genetics
from Bio import AlignIO
import numpy as np
import os
ZERO=1e-8
g = Genetics()


def build_mu_dict(pi, kappa):
    ''' Construct mutation rates dictionary.'''
    a = pi[0]
    c = pi[1]
    g = pi[2]
    t = pi[3]
    mu = {'AG':kappa*g, 'TC':kappa*c, 'GA':kappa*a, 'CT':kappa*t, 'AC':c, 'TG':g, 'CA':a, 'GT':t, 'AT':t, 'TA':a, 'GC':c, 'CG':g}  
    return mu


def extract_mutation_rates(hyphile): # hehe argument!

    with open(hyphile, "r") as f:
        hyout_string = f.read()
        
    # Nucleotide frequencies #
    nuc_freqs = np.zeros(4)
    find_aux_1 = re.search("aux_1=(0\.\d+)", hyout_string)
    aux_1 = float( find_aux_1.group(1) )
    find_aux_2 = re.search("aux_2=(0\.\d+)", hyout_string)
    aux_2 = float( find_aux_2.group(1) )
    find_aux_3 = re.search("aux_3=(0\.\d+)", hyout_string)
    aux_3 = float( find_aux_3.group(1) )
    nuc_freqs[0] = aux_1
    nuc_freqs[1] = (1. - aux_1) * aux_2
    nuc_freqs[2] = (1. - aux_1) * (1. - aux_2) * aux_3
    nuc_freqs[3] = (1. - aux_1) * (1. - aux_2) * (1. - aux_3)
    assert( abs(1. - np.sum(nuc_freqs)) < ZERO ), "Nucleotide frequencies do not sum to 1!"
    nuc_freqs_ordered = nuc_freqs[[3,1,0,2]]

    # Kappa #
    find_kappa = re.search("k=(\d+\.\d+)", hyout_string)
    kappa = find_kappa.group(1)
    
    # Mutation dictionary #
    return build_mu_dict(nuc_freqs, float(kappa) )



def aa_to_codon_freqs(aa_freqs):
    
    codon_freqs = np.zeros(61)
    aa_dict = dict(zip(g.amino_acids, aa_freqs))
    for aa in aa_dict:
        syn_codons = g.genetic_code[ g.amino_acids.index(aa) ]
        cf = aa_dict[aa] / float(len(syn_codons))
        for syn in syn_codons:
            codon_freqs[ g.codons.index(syn) ] = cf
    
    assert(1. - np.sum(codon_freqs) <= ZERO), "bad codon freqs from aa freqs"
    return codon_freqs





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
    
    # Ignore gaps or noncanonical in column
    column = []
    print "\n\n"
    print site
    column_raw = list(records[:, site:site+3])
    for entry in column_raw:
        codon_current = str(entry.seq)
        print codon_current
        assert(codon_current not in g.stop_codons), "WOAH, STOP CODON!"
        if codon_current in g.codons:
            column.append(codon_current)
    length = len(column)
    
    for codon in column:
        freq_dict[codon] += 1./length
    
    for aa_family in g.genetic_code:
        total_aa = 0.
        for entry in aa_family:
            total_aa += freq_dict[entry]
        total_aa /= float(len(aa_family))
        for entry in aa_family:
            freq_dict_from_aa[entry] = total_aa
           
    dnds_c = compute_dnds(freq_dict, mu_dict)       
    dnds_a = compute_dnds(freq_dict_from_aa, mu_dict)       

    return dnds_c, dnds_a

#, 
#stop "PF02518", "PF01336", "PF01926","PF01266",
emp_datasets = ["PF00126", "PF04055", "PF00593","PF07715"] #["lysin","camelid", "pb2" ,"amine"]
sim_datasets = ["1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1RII_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
datasets = {"empirical":emp_datasets} #"simulation":sim_datasets

for datatype in datasets:
    
    seqdir = "../data/" + datatype + "/"
    resdir = "../results/" + datatype + "/"
    
    for data in datasets[datatype]:
        print data
        
        outfile = resdir + data + "_dnds_from_alignment.txt"
        alnfile = seqdir + data + ".fasta"
        hyphile = resdir + data + "_swmutsel_hyout.txt"

        mu_dict = extract_mutation_rates(hyphile)
    
        records = AlignIO.read(alnfile, "fasta")
        omegas = []
    
        for i in range(len(records[0])/3):
            site = i*3
            omegas.append( derive_site_dnds(site, records, mu_dict) )
    
        
        with open(outfile, "w") as outf:
            outf.write("site,dnds_codon,dnds_aa\n")
            for i in range(len(omegas)):
                outf.write(str(i) + "," + str(omegas[i][0]) + "," + str(omegas[i][1]) + "\n")
    
                
                        

