# SJS
# This script computes Jensen-Shannon divergence between inferred frequencies and true frequencies, for simulated datasets 

import os
import numpy as np
from universal_functions import *


def calculate_save_jsd(true, fitness, mu_dict, outfile):
    
    jsd = []
    for i in range(len(fitness)):
        # convert codon frequencies to amino acid frequencies
        truesite = np.array(codon_freqs_to_aa_freqs(true[i]))

        # inferred
        raw_inferred  = np.array( codon_freqs_from_fitness_eigenvector(fitness[i], mu_dict) )
        inferredsite = codon_freqs_to_aa_freqs(raw_inferred)

        assert(len(truesite) == 20), "\nTrue aa frequencies are not length 20?"
        assert(len(inferredsite) == 20), "\nInferred aa frequencies are not length 20?"
       
        jsd.append( calculate_jsd(truesite, inferredsite) )     
        
    assert(len(jsd) == len(true)), "\nJSD not fully calculated for a dataset."
    with open(outfile, "w") as outf:
        outf.write("site,jsd\n")
        outf.write( "\n".join([str(i+1)+","+str(jsd[i]) for i in range(len(jsd))]) )

  
    

def calc_kl(a,b):
    return np.sum( a * np.log(a/b) ) 

    
### Compute JSD from frequency distributions ###
def calculate_jsd(p, q):
    '''
        Calculate JSD between frequency distributions.
        p and q are input distributions to compare, each of length 20.
    '''

    m = (p+q)/2.
    
    term1 = 0.5*calc_kl(p, m)
    term2 = 0.5*calc_kl(q, m)
    
    return np.sqrt( term1 + term2 )




def main():
    
    method_suffixes = {"swmutsel":"_MLE.txt", "phylobayes":".aap"}
    
    rawdir = "../results/raw_results/"
    outdir = rawdir + "derived_dnds_coeffs/"
    
    truedir = "simulation/true_simulation_parameters/"
    
    for method in method_suffixes:

        indir = rawdir + method + "/"
        suffix = method_suffixes[method] 
        
        filenames = os.listdir(indir)
        for file in filenames:

            prefix = None
            if file.endswith(suffix):
                prefix = file.split(suffix)[0]
                simprefix = "_".join( prefix.split("_")[:-1] )
                outprefix = prefix.replace("_simulated", "")
                simprefix = simprefix.replace("_simulated","")
                
                outfile = outdir + outprefix + "_jsd.csv"
    
                inf_fitnesses, mu_dict = extract_parameters(indir, prefix)                     
                true_codon_frequencies = np.loadtxt(truedir + simprefix + "_true_codon_frequencies.txt")
                
                if not os.path.exists(outfile):
                    print "Computing JSD for", outprefix
                    calculate_save_jsd(true_codon_frequencies, inf_fitnesses, mu_dict, outfile)

main()      