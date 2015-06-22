# SJS
# Here, we predict site-specific dN/dS values from the frequencies in the alignment columns (using true kappa)

import os
import sys
sys.path.append("/Users/sjspielman/Research/sitewise_dnds_mutsel_exp/scripts/simulation")
from dnds_mutsel_functions import *
import numpy as np
from Bio import AlignIO
import multiprocessing

# GLOBALS, HUZZAH
numproc = 3 # processors
resdir = "/Users/sjspielman/Research/mutsel_bench/results/"
seqdir = "/Users/sjspielman/Research/mutsel_bench/data/"

mu = 1.
kappa = 4.0
mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]



def compute_dnds(d):
    dnds = None
    # Assertion error will be raised in dnds calculation if no evolution has occurred. These are uninformative sites.
    try:
        dnds = derive_dnds( d, mu_dict )
        if np.isnan(dnds) or np.isinf(dnds):
            dnds = "NA"
    except AssertionError:
        dnds = "NA"
    
    assert(dnds != None), "dnds not computed."
    return dnds


def multi_worker(name):
    '''
        Wrapper function for computing dnds in multiprocessing framework.
    '''
    alnfile = seqdir + name + ".phy"
    outfile = resdir + name + "_predicted_dnds.csv"
    records = AlignIO.read(alnfile, "phylip-relaxed")
    omegas = [] # list of lists, where first entry is dnds from codon frequencies and second entry is dnds from amino acid frequencies
    
    for i in range(0, len(records[0]), 3):
        freq_dict = dict(zip(codons, np.repeat(0., 61)))
        freq_dict_from_aa = {}  # dict(zip(codons, np.repeat(0., 61)))
        
        column = list(records[:, i:i+3])
        length = float(len(column))
        
        for row in column:
            codon = str(row.seq)
            freq_dict[codon] += 1./length
        
        for aa_family in genetic_code:
            total_aa = 0.
            for entry in aa_family:
                total_aa += freq_dict[entry]
            total_aa /= float(len(aa_family))
            for entry in aa_family:
                freq_dict_from_aa[entry] = total_aa
        dnds1 = compute_dnds(freq_dict)
        dnds2 = compute_dnds(freq_dict_from_aa)
        omegas.append( [dnds1, dnds2] )                     
    
    with open(outfile, "w") as outf:
        outf.write("site,dnds_codon,dnds_aa\n")
        for i in range(len(omegas)):
            outf.write(str(i) + "," + str(omegas[i][0]) + "," + str(omegas[i][1]) + "\n")

# 
# roots = ["rootrandom", "rootmax", "rootmin"]
# syns  = ["synsel", "nosynsel"]
# ns    = ["n5","n6","n7","n8", "n9", "n10"]
# bls   = ["bl0.005", "bl0.01", "bl0.02", "bl0.04", "bl0.08", "bl0.16"]
# job_names = []
# total_runs = len(roots) * len(syns) * len(ns) * len(bls)
# 
# 
# 
# for roottype in roots:
#     for s in syns:
#         for n in ns:
#             for bl in bls:
#                 name = n + "_" + bl + "_" + s + "_" + roottype
#                     job_names.append(name)
# 
job_names = ["n9_bl0.16"]
total_runs = 1
completed_runs = 0
while (completed_runs < total_runs):
    if (total_runs - completed_runs) >= numproc:
        jobs = []
        for i in range(completed_runs, completed_runs + numproc):
            name = job_names[i]
            print name
            p = multiprocessing.Process(target = multi_worker, args = (name,))
            jobs.append(p)
            p.start()
        for j in jobs:
            j.join()
        completed_runs += numproc
    
    else:
        jobs = []
        for i in range(total_runs - completed_runs - 1, total_runs):
            name = job_names[i]
            p = multiprocessing.Process(target = multi_worker, args = (name,))
            jobs.append(p)
            p.start()
        for j in jobs:
            j.join()
        break            
         

                
                
                
                
                
                
                
                
                        

