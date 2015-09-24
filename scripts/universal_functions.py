# SJS
# This file contains functions which are generally used throughout analyses.


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




### Compute dN/dS from MutSel parameters ###
def dnds_from_params(site_fitness, mu_dict):
    
    # Build the MutSel matrix (assumes equal codon frequencies per amino acid) with these parameters and extract equilibrium codon frequencies
    # Note that if mutation is symmetric, this gives the same results as Boltzmann distribution does anyways
    eqfreqs = codon_freqs_from_fitness_eigenvector(site_fitness, mu_dict)
    
    # Derive dN/dS and return
    c = dNdS_from_MutSel( dict(zip(g.codons, eqfreqs)), mu_dict)
    return c.compute_dnds()
 

def jsd_from_params(true_codon_freqs, site_fitness, mu_dict):
    
    # convert codon frequencies to amino acid frequencies
    true = np.array(codon_freqs_to_aa_freqs(true_codon_freqs)) + ZERO   # since need to log
    
    # inferred
    raw_inferred  = codon_freqs_from_fitness_eigenvector(site_fitness, mu_dict)
    inferred = np.array(codon_freqs_to_aa_freqs(raw_inferred)) + ZERO   # since need to log

    # JSD between
    return calculate_jsd( true, inferred )
    

def calc_kl(a,b):
    return np.sum( a * np.log(a/b) ) 

    
### Compute JSD from frequency distributions ###
def calculate_jsd(p, q):
    '''
        p and q are input distributions to compare, each of length 20.
    '''

    m = (p+q)/2.


    term1 = 0.5*calc_kl(p, m)
    term2 = 0.5*calc_kl(q, m)
    
    return np.sqrt( term1 + term2 )



########## Generic functions for converting among codon, amino-acid frequencies and fitnesses ###########


def codon_freqs_from_fitness_boltzmann(codon_fitness):
    '''
        Convert codon fitness to equilibrium frequencies using Boltzmann distribution (Sella and Hirsh 2005).
    '''
    codon_freqs = np.zeros(61)
    count = 0
    for codon in g.codons:
        codon_freqs[count] = np.exp( codon_fitness[codon] )
        count += 1
    codon_freqs /= np.sum(codon_freqs)                   
    assert( abs(1. - np.sum(codon_freqs)) <= ZERO), "codon_freq doesn't sum to 1 in codon_fitness_to_freqs"
    return codon_freqs 


def codon_freqs_from_fitness_eigenvector(fitness, mu):
    '''
        Extract equilibrium frequencies from eigenvector of MutSel matrix.
    '''
    params = {"fitness": fitness, "mu": mu}
    m = Model("mutsel", params)
    matrix = m.matrix
    eqfreqs = m.extract_state_freqs()
    return eqfreqs



def aa_freqs_to_codon_freqs(aa_freqs):
    '''
        Convert amino-acid frequencies to codon frequencies, assuming no bias.
    '''
    codon_freqs = np.zeros(61)
    aa_dict = dict(zip(g.amino_acids, aa_freqs))
    for aa in aa_dict:
        syn_codons = g.genetic_code[ g.amino_acids.index(aa) ]
        cf = aa_dict[aa] / float(len(syn_codons))
        for syn in syn_codons:
            codon_freqs[ g.codons.index(syn) ] = cf
    
    assert(1. - np.sum(codon_freqs) <= ZERO), "bad codon freqs from aa freqs"
    return codon_freqs, dict(zip(g.codons, codon_freqs))


def codon_freqs_to_aa_freqs(codonfreqs):
    '''
        Codon codon frequencies to amino-acid frequencies.
    '''
    cf = dict(zip(g.codons, codonfreqs))
    aa_freqs = []
    for family in g.genetic_code:
        total = 0.
        for syn in family:
            total += cf[syn]
        total /= float(len(family))
        aa_freqs.append(total)
        
    return aa_freqs


def aa_fitness_to_codon_fitness(fitness):
    '''
        Convert amino-acid fitnesses to codon fitnesses.
    '''
    
    d = {}
    for i in range(20):
        syn_codons = g.genetic_code[i]
        for syn in syn_codons:
            d[ syn ] = fitness[i]   
    codon_fitness = np.zeros(61)
    count = 0
    for i in range(61):
        codon_fitness[i] = d[g.codons[i]]
    return codon_fitness


#########################################################################################



########## Parsing functions for swMutSel and pbMutSel inferences to get dN/dS, JSD ###########



def build_mu_dict(pi, kappa):
    ''' Construct mutation rates dictionary.'''
    a = pi[0]
    c = pi[1]
    g = pi[2]
    t = pi[3]
    mu = {'AG':kappa*g, 'TC':kappa*c, 'GA':kappa*a, 'CT':kappa*t, 'AC':c, 'TG':g, 'CA':a, 'GT':t, 'AT':t, 'TA':a, 'GC':c, 'CG':g}  
    return mu
    



def parse_swMutSel_mutation(infile):
    '''
        Extract kappa (line 11) and nucleotide frequencies (line 14) from a swMutSel MLE file.
        Combine info to return a dictionary of mutation rates.
    '''
    with open(infile, "r") as f:
        lines = f.readlines()
    kappa = float(lines[10].strip())
    rawpis   = lines[13].strip().split(",")
    pis = []
    for pi in rawpis:
        pis.append(float(pi)) 
    
    return build_mu_dict(pis, kappa)

    


def parse_pbMutSel_fitness(cpu, job_name, burnin = '100'):
    '''
        Parse aap and tracefile. Return site-wise amino acid fitnesses and dictionary of mutation rates.
        Fitness values are return as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ], and return dictionary of mutation rates.
    '''
    # Read in fitness values from .aap file and take exponential to obtain usable fitness values
    fitness_raw = np.loadtxt(job_name + ".aap")
    fitness = np.exp(fitness)    
    
    # Grab posterior means for mutation rates from the .trace file, using same sampling as for fitnesses
    mu_dict = parse_trace_mutation( job_name + ".trace", burnin )
    
    return fitness, mu_dict



def parse_pbMutSel_mutation(file, burnin = 100):
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
    
