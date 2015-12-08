# SJS
# This script computes dN/dS values from MutSel parameters (both sitewise fitness and mutation rates) as inferred with swMutSel and pbMutSel.

import sys
import numpy as np
from universal_functions import *


def extract_info(directory, name):
    '''
        Extract fitness and mutation rates from a given inference.
    '''
    fitness = None
    mu_dict = None
    if "phylobayes" in name:
        fitness = np.loadtxt(directory + name + ".aap")
        fitness = np.log(fitness)
        mu_dict = parse_pbMutSel_mutation(directory + name + ".trace")
    else:
        fitness = np.loadtxt(directory + name + "_fitness.txt")
        mu_dict = parse_swMutSel_mutation(directory + name + "_MLE.txt")   

    assert(fitness is not None), "\n Could not retrieve fitness values."
    assert(mu_dict is not None), "\n Could not retrieve mutation rates."
    return fitness, mu_dict
    
    

def calculate_save_dnds(fitness, mu_dict, outfile):
    '''
        Perform site-wise dN/dS calculations.
    '''
    dnds = []
    for sitefitness in fitness:
        dnds.append( dnds_from_params(sitefitness, mu_dict) )     
    with open(outfile, "w") as outf:
        outf.write( "\n".join([str(i) for i in dnds]) )




def calculate_save_coeffs(fitness, outfile):
    '''
        Compute and save distribution of selection coefficients.
    '''
    all_coefficients = []
    for site in fitness:
        for i in range(len(site)):
            f_i = site[i]
            for j in range(len(site)):
                if i == j:
                    continue
                else:
                    all_coefficients.append(f_i - site[j])
    with open(outfile, "w") as outf:
        outf.write( "\n".join([str(i) for i in all_coefficients]) )
    
    


def compute_dnds_coefficicents(input_directory, output_directory):
    '''
        For all inference files in a given input directory, calculate dN/dS and selection coefficients.
    '''
    suffixes = ["_MLE.txt", ".aap"]
    filenames = os.listdir(input_directory)
    for file in filenames:
    
        prefix = None
        for suffix in suffixes:
            if file.endswith(suffix):
                prefix = file.split(suffix)[0]
            
                outfile_dnds   = output_directory + prefix + "_dnds.txt"
                outfile_coeffs = output_directory + prefix + "_selcoeffs.txt"
                
                fitness, mu_dict = extract_info(input_directory, prefix)
                
                if not os.path.exists(outfile_dnds):
                    print "Computing dN/dS for", prefix
                    calculate_save_dnds(fitness, mu_dict, outfile_dnds)
                
                if not os.path.exists(outfile_coeffs):
                    print "Computing selection coefficients for", prefix
                    calculate_save_coeffs(fitness, outfile_coeffs)


def main():
    
    # Simulated yeast and empirical inferences
    for datatype in ["simulation", "empirical"]:
        rawdir = "../results/raw_results/" + datatype + "/"
        outdir = rawdir + "derived_dnds_coeffs/"
    
        for method in ["phylobayes", "swmutsel"]:
            indir = rawdir + method + "/"
            compute_dnds_coefficicents(indir, outdir)


    # Simulated from empirical 
    rawdir = "../results/raw_results/simulated_from_empirical/"
    indir = rawdir + "swmutsel/"
    outdir = rawdir + "derived_dnds_coeffs/"
    compute_dnds_coefficicents(indir, outdir)


main()

    
    