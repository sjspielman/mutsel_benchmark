# SJS
# This script computes Jensen-Shannon divergence between inferred frequencies and true frequencies, for simulated datasets 

import sys
import numpy as np
from scipy import stats
from universal_functions import *





datasets = ["1B4T_A", "1G58_B", "1GV3_A", "1HUR_A", "1IBS_A", "1PV1_A", "1QMV_A", "1R6M_A", "1V9S_B", "1W7W_B", "1X1O_B", "1YPI_A", "1ZNN_A", "2A84_A", "2BCG_Y", "2CFE_A", "2CJM_C", "2CNV_A", "2FLI_A", "2G0N_B"]
methods  = ["phylobayes"] #, "mvn10", "mvn100", "mvn1000", "d1.0", "d0.1", "d0.01", "nopenal"]
resdir = "../results/simulation/"
truedir = "simulation/flib/"
    
for data in datasets:  
    print data
    true_codon_frequencies = np.loadtxt(truedir + data + "_codon_freq_lib.txt") 

    for method in methods:
            
        outfile = resdir + data + "_" + method + "_jsd.txt"
        print "Calculating JSD for", data, method
            
        ######## OBTAIN FITNESS VALUES AND MUTATION RATES ########
        if method == "phylobayes":
            fitness = np.loadtxt(resdir + data + "_simulated_phylobayes.aap")
            fitness = np.log(fitness)
            mu_dict = parse_pbMutSel_mutation(resdir + data + "_simulated_phylobayes.trace")
        else:
            fitness = np.loadtxt(resdir + data + "_simulated_" + method + "_fitness.txt")
            mu_dict = parse_swMutSel_mutation(resdir + data + "_simulated_" + method + "_MLE.txt")   


        ######## CALCULATE AND SAVE JSD ########
        jsd = []
        count = 0
        for sitefitness in fitness:
            jsd.append( jsd_from_params(true_codon_frequencies[count], sitefitness, mu_dict) ) 
            count += 1    
        with open(outfile, "w") as outf:
            outf.write( "\n".join([str(i) for i in jsd]) )
            







