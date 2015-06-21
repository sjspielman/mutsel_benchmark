# SJS 2/10/15. 
# Parse phylobayes results to get dN/dS.

import subprocess
import sys
import numpy as np
from dnds_functions import *



def parse_pbout( job_name, cpu, burnin = "100", until = "1100" ):
    ''' Create the .aap file to grab amino acid fitnesses from. '''


    readpb_call = "mpirun -np " + str(cpu) + " ./readpb_mpi -x " + str(burnin) + " 1 " + str(until) + " " + job_name + "\n"
    print readpb_call
    run_readpb_call = subprocess.call(readpb_call, shell = True)
    assert( run_readpb_call == 0 ), "readpb_mpi didn't run!"


def parse_trace_mutation( job_name, burnin = "100" ):   
    ''' 
        Parse the .trace file returned by phylobayes to obtain mutation rates
        ORDER (starts at field 14, index from 0): nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT
    '''
    params = np.loadtxt(job_name + ".trace", skiprows = int(burnin) + 1) # skip header + burnin
    params = params[::,14::] # Sample rows according to every and keep only the relevant columns
    post_mean = np.mean(params, axis = 0) # Posterior means for desired parameters
    mu_dict = {'AC': post_mean[0], 'CA': post_mean[0], 'AG': post_mean[1], 'GA': post_mean[1], 'AT': post_mean[2], 'TA': post_mean[2], 'CG': post_mean[3], 'GC': post_mean[3], 'CT': post_mean[4], 'TC': post_mean[4], 'GT': post_mean[5], 'TG': post_mean[5] }
    
    return mu_dict



def main():
    
    
    usage = "\nUsage: python process_phylobayes.py <job_name> <cpu>. Note that all files and all executables ('pb_mpi' and 'readpb_mpi') must be in the working directory!"
    assert( len(sys.argv) == 3 ), usage
    job_name = sys.argv[1]
    cpu = sys.argv[2]
    
 
    # Read in fitness values from .aap file
    parse_pbout( job_name, cpu ) # single cpu for postprocessing.
    sitewise_fitness = np.loadtxt(job_name + ".aap")
    
    # Grab posterior means for mutation rates from the .trace file, using same sampling as for fitnesses
    mu_dict = parse_trace_mutation( job_name, cpu )
    
    # Compute site-wise dN/dS from amino-acid fitness values
    dnds = []
    for site_fitness in sitewise_fitness:
        cf = fitness_to_codonfreq(site_fitness)
        dnds.append( derive_dnds( cf, mu_dict ) )
    
    # Save dnds file
    with open(job_name + "_dnds.txt", "w") as outf:
        outf.write( "\n".join([str(i) for i in dnds]) )
    
    
    
main()
    


