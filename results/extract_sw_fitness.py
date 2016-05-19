import os
from numpy import savetxt

files = [x for x in os.listdir("swmutsel/") if x.endswith("MLE.txt")]

for file in files:
    fitfile = file.split("_MLE.txt")[0] + "_fitness.txt"
    if os.path.exists("swmutsel/" + fitfile):
        continue
    else:
        print file
        fitness = []
        new_order = [0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17] # need to reorder fitnesses from tamuri's mapping to mine.
        with open("swmutsel/" + file, 'r') as f:
            for line in f:
                newline = line.split(',')[1:] # first col in csv is site index, so ignore it.
                if len(newline) == 20:
                    fitness.append( [float(y) for (x,y) in sorted(zip(new_order,newline))] )

        # Save a file with the correctly-ordered (well, my order) amino-acid fitness values
        savetxt("swmutsel/" + fitfile, fitness, delimiter = '\t')        
