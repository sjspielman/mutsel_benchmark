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
        fitness = np.loadtxt(directory + name + "_phylobayes.aap")
        fitness = np.log(fitness)
        mu_dict = parse_pbMutSel_mutation(directory + name + "_phylobayes.trace")
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


def compute_all_dnds(input_directory, output_directory):
    '''
        For all inference files in a given input directory, calculate dN/dS.
    '''
    suffixes = ["_MLE.txt", ".aap"]
    filenames = os.listdir(input_directory)
    for file in filenames:
    
        prefix = None
        for suffix in suffixes:
            if file.endswith(suffix):
                prefix = file.split(suffix)[0]
            
                outfile = output_directory + prefix + "_dnds.txt"
                if os.path.exists(outfile):
                    continue
            
                fitness, mu_dict = extract_info(input_directory, prefix)
                print "Calculating dN/dS for simulated", prefix
                calculate_save_dnds(fitness, mu_dict, outfile)


def main():
    
    # Simulated yeast and empirical inferences
    for datatype in ["simulation", "empirical"]:
        rawdir = "../results/raw_results/" + datatype + "/"
        outdir = rawdir + "derived_dnds/"
    
        for method in ["phylobayes", "swmutsel"]:
            indir = rawdir + method + "/"
            compute_all_dnds(indir, outdir)


    # Simulated from empirical 
    rawdir = "../results/raw_results/simulated_from_empirical/"
    indir = rawdir + "swmutsel/"
    outdir = rawdir + "derived_dnds/"
    compute_all_dnds(indir, outdir)


main()

    
    