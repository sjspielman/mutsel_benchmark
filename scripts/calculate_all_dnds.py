# SJS
# This script computes dN/dS values from MutSel parameters (both sitewise fitness and mutation rates) as inferred with swMutSel and pbMutSel.

import sys
import numpy as np
from Bio import AlignIO
from compute_dnds_from_mutsel import *
from universal_functions import *


    

#"1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", 

emp_datasets = ["amine", "pb2", "PF00593", "PF01266", "PF01336", "PF01926" ,"PF00126", "PF04055", "PF00593","PF07715"]
sim_datasets = ["1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
datasets = {"empirical":emp_datasets, "simulation":sim_datasets}
inftypes = {"empirical": ["mvn10"], "simulation": ["mvn10", "mvn100", "mvn1000", "d1.0", "d0.1", "d0.01", "nopenal"]}

for datatype in datasets:
    
    seqdir = "../data/" + datatype + "/"
    resdir = "../results/" + datatype + "/"
    
    for data in datasets[datatype]:
        print data
        
        for method in inftypes[datatypes]:
            print method
            
            ######## OBTAIN FITNESS VALUES AND MUTATION RATES ########
            if method == "pbmutsel":
                fitness = np.loadtxt(prefix + "_phylobayes.aap")
                tracefile = resdir + data + "_phylobayes.trace"
                mu_dict = parse_pbMutSel_mutation(tracefile)
            else:
                fitness = np.loadtxt(resdir + data + "_" + method + "_fitness.txt")
                mlefile = resdir + data + "_" + method + "_MLE.txt"
                mu_dict = parse_swMutSel_mutation(mlefile)   
            
            

            ######## CALCULATE AND SAVE dN/dS ########
            dnds = []
            for sitefitness in fitness:
                dnds.append( dnds_from_params(sitefitness, mu_dict) )
            outfile = resdir + data + "_" + method + "_dnds.txt"     
            with open(outfile, "w") as outf:
                outf.write( "\n".join([str(i) for i in dnds]) )
            
