# SJS
# Given a coding-sequence alignment and phylogeny, run the PhyloBayes mutsel model (Rodrigue 2010, 2014). 
# Once site-wise fitness values are inferred, compute site-wise dN/dS values.

# Usage: python run_phylobayes_mutsel.py <alignment_file> <tree_file> <cpu> <job_name> . 
## Alignment must be in phylip format, and tree must be in newick format. Both must be in working directory.
## Executables pb_mpi and readpb_mpi must be in the working directory.


import re
import os
import sys
import subprocess
import numpy as np
from Bio import AlignIO
from dendropy import Tree
from pyvolve import *
from compute_dnds_from_mutsel import *
from parsing_functions import *
g = Genetics()
 

def main():
    
    restart = bool(sys.argv[1])
    cpu = sys.argv[2]
    job_name = sys.argv[3]
    
    try:
        alnfile = sys.argv[4]
    except:
        assert(restart is True), "Specified alignment file does not exist. Path?"
    try:
        treefile = sys.argv[5]
    except:
        assert(restart is True), "Specified tree file does not exist. Path?"
    
    
    
    if not restart:
        # Rewrite tree to create trifurcating root, as needed by phylobayes mpi
        tree = Tree.get_from_path(treefile, "newick", rooted = False)
        tree.resolve_polytomies() # in case of polytomies.
        tree.update_splits() # this will create a trifurcating root on an unrooted tree
        tstring = str(tree).replace('[&U] ', '')
        with open('temp.tre', 'w') as tf:
            tf.write(tstring + ';\n')
            
        pb_call = "mpirun -np " + str(cpu) + " ./pb_mpi -mutsel -cat -d " + alnfile + " -T " + treefile + " -x " + str(every) + " " + str(until) + " " + job_name
    
    else:
        pb_call = "mpirun -np 16 " + str(cpu) + " ./pb_mpi " + job_name    
    
    print pb_call
    run_pb_call = subprocess.call(pb_call, shell = True)
    assert( run_pb_call == 0 ), "pb_mpi didn't run!"

    sitewise_fitness, mu_dict = parse_pb(cpu, job_name)

    # Compute site-wise dN/dS from amino-acid fitness values. 
    # For this, we need to build a MutSel model to obtain equilibrium frequences (as mu not symmetric) from left eigenvector
    dnds = []
    for site_fitness in sitewise_fitness:
        dnds.append(dnds_from_params(site_fitness, mu_dict))
    
    
    # Save dnds file
    with open(job_name + "_dnds.txt", "w") as outf:
        outf.write( "\n".join([str(i) for i in dnds]) )
    
    
    
main()
    


