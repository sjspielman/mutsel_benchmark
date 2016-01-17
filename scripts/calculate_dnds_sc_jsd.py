# SJS
# This script computes dN/dS ratios, distributions of selection coefficients, and  from MutSel parameters (both sitewise fitness and mutation rates) as inferred with swMutSel and pbMutSel.

import sys
import numpy as np
from universal_functions import *



  
def calc_kl(a,b):
    '''
        Compute KL distance from b to a.
    '''
    return np.sum( a * np.log(a/b) ) 

    

def calculate_jsd(p, q):
    '''
        Calculate JSD between frequency distributions.
        p and q are input distributions to compare, each of length 20.
    '''

    m = (p+q)/2.
    
    term1 = 0.5*calc_kl(p, m)
    term2 = 0.5*calc_kl(q, m)
    
    return np.sqrt( term1 + term2 )
    
    
def calculate_save_jsd(true, fitness, mu_dict, outfile):
    '''
        Calculate JSD between inferred and true frequency distribution and save.
    '''    
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



def calculate_save_dnds(fitness, mu_dict, outfile):
    '''
        Perform site-wise dN/dS calculations and save.
    '''
    dnds = []
    for sitefitness in fitness:
        dnds.append( dnds_from_params(sitefitness, mu_dict) )     
    with open(outfile, "w") as outf:
        outf.write("site,dnds\n")
        outf.write( "\n".join([str(i+1)+","+str(dnds[i]) for i in range(len(dnds))]) )
    
    return dnds




def calculate_save_coeffs(fitness, outfile, dnds):
    '''
        Compute and save distribution of selection coefficients.
    '''
    raw = []
    binned = []
    sitednds = []
    y = 0
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
                    sitednds.append(dnds[y])
        y += 1
    with open(outfile, "w") as outf:
        outf.write("realcoeff,binnedcoeff,sitednds\n")
        outf.write( "\n".join([str(raw[x])+","+str(binned[x])+","+str(sitednds[x]) for x in range(len(raw))]) )
    
    




def compute_values(input_directory, output_directory, true_directory, suffix):
    '''
        For all inference files in a given input directory, calculate dN/dS, JSD, and selection coefficients.
    '''
    
    filenames = os.listdir(input_directory)
    for file in filenames:
    
        prefix = None
        if file.endswith(suffix):
            prefix = file.split(suffix)[0]
            
            outfile_dnds   = output_directory + prefix + "_dnds.csv"
            outfile_coeffs = output_directory + prefix + "_selcoeffs.csv"
            outfile_jsd    = output_directory + prefix + "_jsd.csv"

            fitness, mu_dict = extract_parameters(input_directory, prefix)                         
            
            # dN/dS and coefficients calculated together, in same if block, because coefficient distributions are saved w/ dnds indicated per-site.
            if not os.path.exists(outfile_dnds) and not os.path.exists(outfile_coeffs):
                print "Computing dN/dS and selection coefficients for", prefix
                dnds = calculate_save_dnds(fitness, mu_dict, outfile_dnds)
                calculate_save_coeffs(fitness, outfile_coeffs, dnds)        
            
            
            if not os.path.exists(outfile_jsd):
                print "Computing JSD for", prefix
                simprefix = "_".join(prefix.split("_")[0:3])
                true_codon_frequencies = np.loadtxt(true_directory + simprefix + "_true_codon_frequencies.txt")
                calculate_save_jsd(true_codon_frequencies, fitness, mu_dict, outfile_jsd)

                





def main():
    
    method_suffixes = {"swmutsel":"_MLE.txt", "phylobayes":".aap"}
    
    truedir = "simulation/true_simulation_parameters/"
    rawdir  = "../results/raw_results/"
    outdir  = rawdir + "dnds_coeffs_jsd/"
    
    for method in method_suffixes:
        indir = rawdir + method + "/"
        suffix = method_suffixes[method] 
        compute_values(indir, outdir, truedir, suffix)

main()

    
    
