# Quick script to extract fitness values from swmutsel output, specifically for simulated_from_empirical results.

import os
import sys
from numpy import savetxt

new_order = [0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17] # need to reorder fitnesses from tamuri's mapping to mine.


resdir = "../results/raw_results/simulation/swmutsel/"
mle_files = os.listdir(resdir)

for file in mle_files:
    if file.endswith("MLE.txt"):
        
        print file
        
        fitness = []
        name = file.split("_MLE.txt")[0]
        fitfile = resdir + name + "_fitness.txt"
        with open(resdir + file, "r") as f:
            for line in f:
                newline = line.split(',')[1:] # first col in csv is site index, so ignore it.
                if len(newline) == 20:
                    fitness.append( [float(y) for (x,y) in sorted(zip(new_order,newline))] )

        savetxt(fitfile, fitness, delimiter = '\t')        



