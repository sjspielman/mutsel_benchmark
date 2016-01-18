# SJS
# Quick script to extract fitness values from inferences. Fitnesses are re-scaled according to true simulated fitness.


import sys
import numpy as np
sys.path.append("../scripts/")
from universal_functions import *


dataset = "1IBS_A"
for deltype in ["strong", "weak"]:
    name = dataset + "_del" + deltype
    outfile = name + "_rescaled_fitness_across_methods.txt" 
    
    inference_directory = "../results/raw_results/"
    true_directory      = "../scripts/simulation/true_simulation_parameters/"
    
    
    # swMutSel
    swmethods = ["nopenal", "mvn1", "mvn10", "mvn100", "d0.01", "d0.1", "d1.0"]

    # Obtain true fitness ranking
    swfit = np.loadtxt(true_directory + name + "_true_aa_fitness.txt")
    site_rescale = [] # The amino acid for this index should have a fitness of 0

    with open(outfile, "w") as f:
        f.write("site,fitindex,fitness,method\n")

    sitecount = 1
    for row in swfit:
        new =  np.copy(row)
        maxfit = np.max(new)
        maxfit_index = int( np.where(new == maxfit)[0][0] )
        site_rescale.append(maxfit_index)
        new -= maxfit
        with open(outfile, "a") as f:
            f.write("\n".join( [str(sitecount) + "," + str(x) + "," + str(new[x]) +",true" for x in range(20)])) 
            f.write("\n")
        sitecount +=1


    for m in swmethods:
        swfit = np.loadtxt(inference_directory + "swmutsel/" + name + "_" + m + "_fitness.txt")
    
        sitecount = 1
        for row in swfit:
            new =  np.copy(row)
            makemax = site_rescale[sitecount - 1]
            new -= new[makemax]
            with open(outfile, "a") as f:
                f.write("\n".join( [str(sitecount) + "," + str(x) + "," + str(new[x]) + "," + m for x in range(20)])) 
                f.write("\n")
            sitecount +=1
        

    pbfit = np.log( np.loadtxt(inference_directory + "phylobayes/" + name + "_phylobayes.aap"))

    sitecount = 1
    for row in swfit:
        new =  np.copy(row)
        makemax = site_rescale[sitecount - 1]
        new -= new[makemax]
        with open(outfile, "a") as f:
            f.write("\n".join( [str(sitecount) + "," + str(x) + "," + str(new[x]) + ",phylobayes" for x in range(20)])) 
            f.write("\n")
        sitecount +=1



