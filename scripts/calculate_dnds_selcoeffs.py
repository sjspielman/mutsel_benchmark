# SJS
# This script computes dN/dS values from MutSel parameters (both sitewise fitness and mutation rates) as inferred with swMutSel and pbMutSel.

import sys
import numpy as np
from universal_functions import *


  

def calculate_save_dnds(fitness, mu_dict, outfile):
    '''
        Perform site-wise dN/dS calculations.
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
    
    


def compute_dnds_coefficicents(input_directory, output_directory, suffix):
    '''
        For all inference files in a given input directory, calculate dN/dS and selection coefficients.
    '''
    
    filenames = os.listdir(input_directory)
    for file in filenames:
    
        prefix = None
        if file.endswith(suffix):
            prefix = file.split(suffix)[0]
            outprefix = prefix.replace("_simulated", "")
            
            outfile_dnds   = output_directory + outprefix + "_dnds.csv"
            outfile_coeffs = output_directory + outprefix + "_selcoeffs.csv"

            fitness, mu_dict = extract_parameters(input_directory, prefix)                     
            
            if not os.path.exists(outfile_dnds) and not os.path.exists(outfile_coeffs):
                print "Computing dN/dS and selection coefficients for", prefix
                dnds = calculate_save_dnds(fitness, mu_dict, outfile_dnds)
                calculate_save_coeffs(fitness, outfile_coeffs, dnds)                
                


def main():
    
    method_suffixes = {"swmutsel":"_MLE.txt", "phylobayes":".aap"}
    
    for datatype in ["simulation", "empirical"]:
        rawdir = "../results/raw_results/" + datatype + "/"
        outdir = rawdir + "derived_dnds_coeffs/"
    
        for method in method_suffixes:
            indir = rawdir + method + "/"
            suffix = method_suffixes[method] 
            compute_dnds_coefficicents(indir, outdir, suffix)

main()

    
    