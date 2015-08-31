# SJS
# Given a coding-sequence alignment and phylogeny, run the swmutsel mutation-selection model (Tamuri 2012,2014) across a variety of penalty functions.
# Can either use HyPhy to optimize mutational parameters and branch lengths before sending to swmutsel, or just optimize branch lengths according to fixed mutation rates.
# Note that optimization is done first with HKY85 to get kappa and nucleotide frequencies, and then with a global MutSel neutral model (this is essentially MG94 with w=1), as this is how swmutsel interprets BL.
# Once site-wise fitness values are inferred, compute site-wise dN/dS values.

# Usage: python run_mutsel.py <alignment_file> <tree_file> <cpu> <dataset>. 
## Alignment must be in either phylip/fasta format, and tree must be in newick format. Both must be in working directory.
## Additionally note that all executables called (swmutsel.jar and HYPHYMP) must be in the working directory.


import re
import os
import sys
import subprocess
import numpy as np
from Bio import AlignIO
from compute_dnds_from_mutsel import *
from pyvolve import *
from parsing_functions import *
g = Genetics()

ZERO=1e-8

# Penalty functions (no penalty + 3 dirichlet + 3 multivariate normal) to run swmutsel with
penalty = {"nopenal" :"", 
           "d0.01"   :" -p dirichlet,0.01", 
           "d0.1"    :" -p dirichlet,0.1", 
           "d1.0"    :" -p dirichlet,1.0", 
           "mvn10"   :" -p mvn,10", 
           "mvn100"  :" -p mvn,100", 
           "mvn1000" :" -p mvn,1000"
          }
    

def run_swmutsel(pi, kappa, alnfile, treefile, cpu, dataset, penalname, penalarg):
    '''
        Call swmutsel and parse output. Return site-wise amino acid fitnesses as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ]
        By default, uses the Dirichlet penalty function with shape 0.1.
        NOTE: Assumes that swmutsel executable is in the current directory and named "swmutsel.jar"   
    '''
    job_name = dataset + '_' + penalname
        
    runstring = "java -jar swmutsel.jar -n " + job_name + " -s " + alnfile + " -t " + treefile + " -fix mutation,branches -k " + kappa + " -pi " + pi + " -gc standard "  + " -threads " + str(cpu) + penalarg
    run_swmutsel = subprocess.call(runstring, shell=True)
    assert(run_swmutsel == 0), "swmutsel did not run!"
    
    fitness = []
    new_order = [0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17] # need to reorder fitnesses from tamuri's mapping to mine.
    with open(job_name + '_MLE.txt', 'r') as f:
        for line in f:
            newline = line.split(',')[1:] # first col in csv is site index, so ignore it.
            if len(newline) == 20:
                fitness.append( [float(y) for (x,y) in sorted(zip(new_order,newline))] )
    
    # Save a file with fitness values
    np.savetxt(job_name + '_fitness.txt', fitness, delimiter = '\t')        
    
    return fitness    









def main():
    usage = "\nUsage: python run_swmutsel.py <aln> <treefile> <cpu>. Note that all files and all executables ('swmutsel'+'HYPHYMP') must be in the working directory!"
    assert( len(sys.argv) == 4 ), usage
    
    dataset = sys.argv[1]
    treefile = sys.argv[2]
    cpu = sys.argv[3]

    alnfile_fasta = dataset + ".fasta"
    alnfile_phy = dataset  + ".phy"
    opt_treefile = dataset + "_swmutsel_optimized.tre"
    
    
    # Prep hyphy input file and call hyphy to optimize mutational parameters, create mu_dict to use later for dnds derivation, and make a treefile with the optimized tree
    subprocess.call("cat " + alnfile_fasta + " " + treefile + " > hyin.txt", shell=True)
    run_hyphy = subprocess.call("./HYPHYMP CPU=" + str(cpu) + " optimize_fmutsel_neutral.bf > hyout.txt", shell=True)
    assert(run_hyphy == 0), "optimize_fmutsel_neutral.bf did not run!"
    
    pi, kappa, treestring, mu_dict = extract_optimized_params("hyout.txt")       
    pi2 = [str(i) for i in pi]
    pi_string = ",".join(pi2)             
    with open(opt_treefile, "w") as f:
        f.write(treestring)

    # Call swmutsel to obtain site-wise fitness values across a variety of penalizations (no penalty and those tested in Tamuri 2014).
    for penal in penalty:
        
        sitewise_fitness = run_swmutsel( pi_string, str(kappa), alnfile_phy, opt_treefile, cpu, dataset, penal, penalty[penal])
    
        # Compute site-wise dN/dS from amino-acid fitness values
        dnds = []
        for site_fitness in sitewise_fitness:
            dnds.append(dnds_from_params(site_fitness, mu_dict))
            
        # Save
        with open(dataset + "_" + penal + "_dnds.txt", "w") as outf:
            outf.write( "\n".join([str(i) for i in dnds]) )
  
    
main()
    


