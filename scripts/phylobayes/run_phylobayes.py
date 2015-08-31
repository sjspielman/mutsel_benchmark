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
g = Genetics()
        

def run_pb(alnfile, treefile, cpu, job_name, every = '5', until = '1100', burnin = '100', restart=False):
    '''
        NOTE: Assumes that pb_mpi and readpb_mpi executables are in the current directory.
    '''
    # Run phylobayes for 5500 generations, sampling every 5, producing posterior distribution of n = 1100
    # use this command if restarting from a trace file
    if restart:
        
    # use this command if from scratch
    else:
        
    


def parse_pb(cpu, job_name, burnin = '100'):
    '''
        Parse tracefile. Return site-wise amino acid fitnesses as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ], and return dictionary of mutation rates.
        Phylobayes chain will have had a length of 5500 and been sampled every 5, resulting in a sample size of 1100. 
        Results are processed here with a burnin of 100, leaving a posterior sample size of 1000. 

    '''
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
    
        # Build the MutSel matrix (assumes equal codon frequencies per amino acid) with these parameters and extract equilibrium codon frequencies
        params = {"mu": mu_dict, "fitness": site_fitness}
        matbuilder = mutSel_Matrix(params, "neutral")
        matrix = matbuilder()
        eqfreqs = matbuilder.extract_state_freqs(matrix)
        cf = dict(zip(g.codons, eqfreqs))
        
        # Derive dN/dS and append to list
        c = dNdS_from_MutSel(cf, mu_dict)
        dnds.append( c.compute_dnds()  )   
    
    
    # Save dnds file
    with open(job_name + "_dnds.txt", "w") as outf:
        outf.write( "\n".join([str(i) for i in dnds]) )
    
    
    
main()
    


