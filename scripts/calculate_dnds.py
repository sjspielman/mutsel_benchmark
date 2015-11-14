# SJS
# This script computes dN/dS values from MutSel parameters (both sitewise fitness and mutation rates) as inferred with swMutSel and pbMutSel.

import sys
import numpy as np
from compute_dnds_from_mutsel import *
from universal_functions import *


emp_datasets =  ["PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715", "pb2"] 
sim_datasets = ["1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
datasets =  {"empirical":emp_datasets, "simulation":sim_datasets}
methods = [{"phylobayes":"phylobayes"}, {"swmutsel": "mvn10", "mvn100", "mvn1000", "d1.0", "d0.1", "d0.01", "nopenal"}]



for datatype in datasets:
    
    resdir_raw = "../results/raw_results/" + datatype + "/"
    outdir = resdir_raw + "derived_dnds/"
    
    for data in datasets[datatype]:        
        for family in methods:

            resdir = resdir_raw + family + "/"   
                  
            for method in family:
                
                outfile = outdir + data + "_" + method + "_dnds.txt"
                print "Calculating dN/dS for", data, method
            
                ######## OBTAIN FITNESS VALUES AND MUTATION RATES ########
                if method == "phylobayes":
                    fitness = np.loadtxt(resdir + data + "_phylobayes.aap")
                    fitness = np.log(fitness)
                    mu_dict = parse_pbMutSel_mutation(resdir + data + "_phylobayes.trace")
                else:
                    fitness = np.loadtxt(resdir + data + "_" + method + "_fitness.txt")
                    mu_dict = parse_swMutSel_mutation(resdir + data + "_" + method + "_MLE.txt")   


                ######## CALCULATE AND SAVE dN/dS ########
                dnds = []
                for sitefitness in fitness:
                    dnds.append( dnds_from_params(sitefitness, mu_dict) )     
                with open(outfile, "w") as outf:
                    outf.write( "\n".join([str(i) for i in dnds]) )
            
