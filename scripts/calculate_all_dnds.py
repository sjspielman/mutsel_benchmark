# SJS
# This script computes dN/dS values from MutSel parameters calulated by swMutSel and pbMutSel.

import sys
import numpy as np
from Bio import AlignIO
from compute_dnds_from_mutsel import *
from universal_functions import *


    



emp_datasets = [] #["PF00126", "PF04055", "PF00593","PF07715"]
sim_datasets = ["1B4T_A_simulated"] #, "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1RII_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
datasets = {"empirical":emp_datasets, "simulation":sim_datasets}
inftypes = ["mvn10", "mvn100", "mvn1000", "d1.0", "d0.1", "d0.01", "nopenal"]#, "pbmutsel"]

for datatype in datasets:
    
    seqdir = "../data/" + datatype + "/"
    resdir = "../results/" + datatype + "/"
    
    for data in datasets[datatype]:
        print data
        
        for method in inftypes:
            print method
            
            ######## OBTAIN FITNESS VALUES ########
            if method == "pbmutsel":
                fitness = np.loadtxt(prefix + "_phylobayes.aap")
            else:
                fitness = np.loadtxt(resdir + data + "_" + method + "_fitness.txt")
            
            
            ######## OBTAIN MUTATION RATES ########
            
            # These inferences have fixed equal mutation rates
            if datatype == "simulation" and method != "pbmutsel":
                sym = True
                mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
            
            # Otherwise, obtain a mutation dictionary from file 
            else:
                sym = False
                mu_dict = obtain_mutation_dictonary(dataset, resdir, method)
            


            ######## CALCULATE AND SAVE dN/dS ########
            dnds = []
            for sitefitness in fitness:
                dnds.append( dnds_from_params(sitefitness, mu_dict) )
            outfile = resdir + data + "_" + method + "_dnds.txt"     
            with open(outfile, "w") as outf:
                outf.write( "\n".join([str(i) for i in dnds]) )
            
