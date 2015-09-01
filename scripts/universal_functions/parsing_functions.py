# Various functions used to parse swMutSel, pbMutSel


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
ZERO=1e-10

        
        

def codon_freqs_from_fitness(fitness, mu):
    params = {"fitness": fitness, "mu": mu}
    matbuilder = mutSel_Matrix(params, scale_matrix = "neutral")
    matrix = matbuilder()
    eqfreqs = matbuilder.extract_state_freqs(matrix)
    return eqfreqs



def dnds_from_params(site_fitness, mu_dict):

    # Build the MutSel matrix (assumes equal codon frequencies per amino acid) with these parameters and extract equilibrium codon frequencies
    eqfreqs = codon_freqs_from_fitness(site_fitness, mu_dict)
    cf = dict(zip(g.codons, eqfreqs))

    # Derive dN/dS and return
    c = dNdS_from_MutSel(cf, mu_dict)
    return c.compute_dnds()
        




def build_mu_dict(pi, kappa):
    ''' Construct mutation rates dictionary.'''
    a = pi[0]
    c = pi[1]
    g = pi[2]
    t = pi[3]
    mu = {'AG':kappa*g, 'TC':kappa*c, 'GA':kappa*a, 'CT':kappa*t, 'AC':c, 'TG':g, 'CA':a, 'GT':t, 'AT':t, 'TA':a, 'GC':c, 'CG':g}  
    return mu
    
    

def extract_optimized_params(infile, return_tree = True):
    '''
        Parse results from hyphy batchfile "optimize_fmutsel_neutral.bf" run.  Return all parameters (a new newick tree string, freqs, kappa). Note that freqs are returned as an ordered (T,C,A,G) csv string and kappa is also a string.
        NOTE: Assumes that the executable HYPHYMP is in the current directory.
    '''    
    with open(infile, "r") as hyout:
        hyout_string = str( hyout.read() )

    # Nucleotide frequencies #
    nuc_freqs = np.zeros(4)
    find_aux_1 = re.search("aux_1=(0\.\d+)", hyout_string)
    aux_1 = float( find_aux_1.group(1) )
    find_aux_2 = re.search("aux_2=(0\.\d+)", hyout_string)
    aux_2 = float( find_aux_2.group(1) )
    find_aux_3 = re.search("aux_3=(0\.\d+)", hyout_string)
    aux_3 = float( find_aux_3.group(1) )
    nuc_freqs[0] = aux_1
    nuc_freqs[1] = (1. - aux_1) * aux_2
    nuc_freqs[2] = (1. - aux_1) * (1. - aux_2) * aux_3
    nuc_freqs[3] = (1. - aux_1) * (1. - aux_2) * (1. - aux_3)
    assert( abs(1. - np.sum(nuc_freqs)) < ZERO ), "Nucleotide frequencies do not sum to 1!"
    nuc_freqs_ordered = nuc_freqs[[3,1,0,2]]

    # Kappa #
    find_kappa = re.search("k=(\d+\.\d+)", hyout_string)
    kappa = find_kappa.group(1)
    
    # Mutation dictionary #
    mu_dict = build_mu_dict(nuc_freqs, float(kappa) )
    
    if return_tree:
        # Tree w/ optimized branch lengths #
        find_tree = re.search("Tree codon_tree=(.+)$", hyout_string)
        tree = find_tree.group(1)
        tree = re.sub(r"Node\d+", r"", tree) # remove node labels
        #tree = re.sub(r"([XN]P)(\d+)(\d)([XN]M)(\d+)(\d)", r"\1_\2.\3_\4_\5.\6", tree) # add back underscores, periods to taxon labels. This is necessary for the amine dataset and shouldn't match anything else.

        return nuc_freqs_ordered, kappa, tree, mu_dict
    
    else:
        return nuc_freqs_ordered, kappa, mu_dict


    


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
    mu_dict = parse_trace_mutation( job_name + ".trace", burnin )
    
    return fitness, mu_dict



def parse_trace_mutation(file, job_name, burnin = 100):
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
    