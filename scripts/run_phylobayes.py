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
from dnds_functions import *

        

def run_pb(alnfile, treefile, cpu, job_name, every = '5', until = '1100', burnin = '100'):
    '''
        Call phylobayes and parse output. Return site-wise amino acid fitnesses as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ]
        By default, phylobayes chain will have a length of 5500, be sampled every 5 (as in Rodrigue 2013 Genetics paper), resulting in a sample size of 1100. Results are then processed with a burnin of 100, leaving a posterior sample size of 1000. 
        NOTE: Assumes that pb_mpi and readpb_mpi executables are in the current directory.
    '''
   
    # Run phylobayes for 5500 generations, sampling every 5, producing n = 1100
    pb_call = "mpirun -np " + str(cpu) + " ./pb_mpi -mutsel -cat -d " + alnfile + " -T " + treefile + " -x " + str(every) + " " + str(until) + " " + job_name
    print pb_call
    run_pb_call = subprocess.call(pb_call, shell = True)
    assert( run_pb_call == 0 ), "pb_mpi didn't run!"
    
    # Parse phylobayes output.
    readpb_call = "mpirun -np " + str(cpu) + " ./readpb_mpi -x " + str(burnin) + " 1 -1 " + job_name + "\n"
    print readpb_call
    run_readpb_call = subprocess.call(readpb_call, shell = True)
    assert( run_readpb_call == 0 ), "readpb_mpi didn't run!"

    # Read in fitness values from .aap file
    fitness = np.loadtxt(job_name + ".aap")
    
    # Grab posterior means for mutation rates from the .trace file, using same sampling as for fitnesses
    mu_dict = parse_trace_mutation( job_name, burnin )
    
    return fitness, mu_dict



def parse_phylobayes_tracefile(job_name, burnin = 100):
    ''' 
        Parse the .trace file returned by phylobayes to obtain mutation rates
        ORDER (starts at field 10, index from 0): fA fC fG fT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT
    '''
   
    # Load trace file, but skip header + burnin and keep only relevant columns (mutational parameters)
    params = np.loadtxt(file, skiprows = int(burnin) + 1, usecols = range(10,20)) 
    post_mean = np.mean(params, axis = 0) # Posterior means for desired parameters
    # Assign post_mean values w/ dictionaries to avoid bugs
    raw_keys = ["A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT"]
    raw_mu   = dict(zip(raw_keys, post_mean)) 
    keys = ["AC", "CA", "AG", "GA", "AT", "TA", "CG", "GC", "CT", "TC", "GT", "TG"]
    mu = {}
    for key in keys:
        targetfreq = raw_mu[ key[1] ]
        exch = raw_mu[ "".join(sorted(key)) ]
        mu[key] = targetfreq * exch
    return mu
    


def main():
    usage = "\nUsage: python run_phylobayes.py <alignment_file> <tree_file> <cpu> <job_name> . Alignment must be in phylip format, and tree must be in newick format. Note that all files and all executables ('pb_mpi' and 'readpb_mpi') must be in the working directory!"
    assert( len(sys.argv) == 5 ), usage
    
    alnfile = sys.argv[1]
    assert( os.path.exists(alnfile) ), "Specified alignment file does not exist. Path?"
    treefile = sys.argv[2]
    assert( os.path.exists(treefile) ), "Specified tree file does not exist. Path?"
    cpu = sys.argv[3]
    job_name = sys.argv[4]
    

    # Rewrite tree to create trifurcating root, as needed by phylobayes mpi
    tree = Tree.get_from_path(treefile, "newick", rooted = False)
    tree.resolve_polytomies() # in case of polytomies.
    tree.update_splits() # this will create a trifurcating root on an unrooted tree
    tstring = str(tree).replace('[&U] ', '')
    with open('temp.tre', 'w') as tf:
        tf.write(tstring + ';\n')
    
    # Call phylobayes to obtain site-wise fitness values and mutation rates
    sitewise_fitness, mu_dict = run_pb( alnfile, 'temp.tre', cpu, job_name )

    # Compute site-wise dN/dS from amino-acid fitness values. 
    # For this, we need to build a MutSel model to obtain equilibrium frequences (as mu not symmetric) from left eigenvector
    dnds = []
    for site_fitness in sitewise_fitness:
    
        # Build the MutSel model with these parameters and extract equilibrium codon frequencies
        matrix = build_mutsel_matrix(mu, codon_fitness)
        eqfreqs = get_eq_from_eig(matrix)
        cf = dict(zip(codons, eqfreqs))
        
        # Derive dN/dS and append to list
        dnds.append( derive_dnds(cf, mu)  )
    
    
    # Save dnds file
    with open(job_name + "_dnds.txt", "w") as outf:
        outf.write( "\n".join([str(i) for i in dnds]) )
    
    
    
main()
    


