# SJS
# Fit HKY85 to each empirical dataset, in order to obtain mutation rates.

import subprocess
import os, shutil

datadir = "../../../data/empirical/"
outdir  = "../../../results/empirical/"

datasets = ["pb2", "PF00126", "PF02518", "PF00593", "PF04055", "PF01266", "PF07715", "PF01336", "PF01926"]

for data in datasets:
    
    print "Fitting HKY85 on", data
    
    # File names
    alnfile = datadir + data + ".fasta"
    treefile = datadir + data + ".tre"
    outfile = outdir + data + "_hky85.fit"
    
    # Set up file for hyphy input
    shutil.copy(alnfile, "hyin.txt")
    os.system("cat " + treefile + " >> hyin.txt")
    
    # Run hyphy
    os.system("HYPHYMP fit_hky85.bf > hyout.txt")
    
    # Save output file
    shutil.move("hyout.txt", outfile)
    
    # Cleanup
    os.system("rm messages.log")
    os.system("rm hyin.txt")
    try:
        os.system("rm errors.log")
    except:
        pass

    
    