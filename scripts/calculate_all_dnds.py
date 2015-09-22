# SJS
# This script computes dN/dS values from MutSel parameters (both sitewise fitness and mutation rates) as inferred with swMutSel and pbMutSel.

import sys
import numpy as np
from Bio import AlignIO
from compute_dnds_from_mutsel import *
from universal_functions import *




emp_datasets = ["amine", "pb2", "PF00126", "PF00593", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715"]
sim_datasets = ["1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
datasets =  {"simulation":sim_datasets, "empirical":emp_datasets}
inftypes = {"empirical": ["mvn10", "phylobayes"], "simulation": ["phylobayes", "mvn10", "mvn100", "mvn1000", "d1.0", "d0.1", "d0.01", "nopenal"]}


for datatype in datasets:
    
    seqdir = "../data/" + datatype + "/"
    resdir = "../results/" + datatype + "/"
    
    for data in datasets[datatype]:        
        for method in inftypes[datatype]:
            
            outfile = resdir + data + "_" + method + "_dnds.txt"
            if os.path.exists(outfile):
                continue
            print "Calculating dN/dS for", data, method
            
            ######## OBTAIN FITNESS VALUES AND MUTATION RATES ########
            if method == "phylobayes":
                fitness = np.loadtxt(resdir + data + "_phylobayes.aap")
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
            
