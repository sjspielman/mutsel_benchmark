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
from copy import deepcopy
from scipy import linalg
from random import uniform, shuffle

g = Genetics()
ZERO=1e-10


### Simulation constants and functions ###
mu_dict        = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
STRONG_FREQ    = 1e-9
WEAK_FIT       = [-6.0, -4.5] # Lowest fitness, highest fitness   
WEAK_THRESHOLD = -10 # Any fitness below here gets increased to an interval in WEAK_FIT
def calculate_save_coeffs(fitness, outfile):
    '''
        Compute and save distribution of selection coefficients.
    '''
    raw = []
    binned = []
    for site in fitness:
        for i in range(len(site)):
            f_i = site[i]
            for j in range(len(site)):
                if i == j:
                    continue
                else:
                    s = f_i - site[j]
                    raw.append(s)
                    if s <= -10:
                        s = -10.
                    if s >= 10.:
                        s = 10.
                    binned.append(s)
    with open(outfile, "w") as outf:
        outf.write("realcoeff,binnedcoeff\n")
        outf.write( "\n".join([str(raw[x])+","+str(binned[x]) for x in range(len(raw))]) )


def apply_weakdel(fitness):
    '''
        Convert a fitness distribution to weakly deleterious regime.
    ''' 
    random_fit = np.random.uniform(low = WEAK_FIT[0], high = WEAK_FIT[1], size = np.sum(fitness <= WEAK_THRESHOLD))
    fitness[fitness <= WEAK_THRESHOLD] = random_fit  
    return fitness



def save_simulation_info(name, frequencies, fitnesses, omegas, entropies):
    '''
        Save simulation parameters to files.
    '''
    freqfile = name + "_true_codon_frequencies.txt"
    fitfile  = name + "_true_aa_fitness.txt"
    selcfile = name + "_true_selcoeffs.csv"
    dnds_entropy_file = name + "_true_dnds_entropy.csv"

    np.savetxt(freqfile, frequencies)
    np.savetxt(fitfile, fitnesses)
    calculate_save_coeffs(fitnesses, selcfile)

    with open(dnds_entropy_file, "w") as f:
        f.write("site,dnds,entropy\n")
        for i in range(len(omegas)):
            f.write(str(i+1)+","+str(omegas[i])+","+str(entropies[i])+"\n")



def calc_entropy(a):
    '''
        Compute entropy of amino acid frequency distribution a.
    '''
    return -1. * np.sum ( a[a > ZERO] * np.log(a[a > ZERO]) )  
    
    
### Compute dN/dS from MutSel parameters ###
def extract_parameters(directory, name):
    '''
        Extract fitness and mutation rates from a given inference.
    '''
    fitness = None
    mu_dict = None
    if "phylobayes" in name:
        fitness = np.loadtxt(directory + name + ".aap")
        fitness = np.log(fitness)
        mu_dict = parse_pbMutSel_mutation(directory + name + ".trace")
    else:
        fitness = np.loadtxt(directory + name + "_fitness.txt")
        mu_dict = parse_swMutSel_mutation(directory + name + "_MLE.txt")   

    assert(fitness is not None), "\n Could not retrieve fitness values."
    assert(mu_dict is not None), "\n Could not retrieve mutation rates."
    return fitness, mu_dict






def dnds_from_params(site_fitness, mu_dict):
    
    # Build the MutSel matrix (assumes equal codon frequencies per amino acid) with these parameters and extract equilibrium codon frequencies
    # Note that if mutation is symmetric, this gives the same results as Boltzmann distribution does anyways
    eqfreqs = codon_freqs_from_fitness_eigenvector(site_fitness, mu_dict)
    
    # Derive dN/dS and return
    c = dNdS_from_MutSel( dict(zip(g.codons, eqfreqs)), mu_dict)
    return c.compute_dnds()
 

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



########## Parsing functions for swMutSel and pbMutSel inferences, specifically for mutation rates ###########
def parse_swMutSel_mutation(infile, return_scaling = False):
    '''
        Extract kappa (line 11) and nucleotide frequencies (line 14) from a swMutSel MLE file.
        Combine info to return a dictionary of mutation rates.
        Argument "return_scaling" means to also return the swMutSel-inferred branch scaling factor.
    '''
    with open(infile, "r") as f:
        lines = f.readlines()
    kappa = float(lines[10].strip())
    rawpis   = lines[13].strip().split(",")
    pis = []
    for pi in rawpis:
        pis.append(float(pi)) 
    t = pis[0]
    c = pis[1]
    a = pis[2]
    g = pis[3]
    mu = {'AG':kappa*g, 'TC':kappa*c, 'GA':kappa*a, 'CT':kappa*t, 'AC':c, 'TG':g, 'CA':a, 'GT':t, 'AT':t, 'TA':a, 'GC':c, 'CG':g}  
    if return_scaling:
        scaling = float(lines[16].strip())
        return mu, scaling
    else:
        return mu

    


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
    
    
def compute_asym(mudict):
    ''' 
        Compute average mu_xy / mu_yx 
    '''
    x = 0  
    completed = []
    ratios = np.zeros(6) 
    for source in ["A", "C", "G", "T"]:
        for target in ["A", "C", "G", "T"]:
            if source == target or source+target in completed or target+source in completed:
                continue
            else:
                ratio = mudict[source + target] / mudict[target + source] 
                if ratio >= 1.:
                    completed.append(source + target)
                    ratios[x] = ratio
                else:
                    completed.append(target + source)
                    ratios[x] = 1./ratio
                x+=1  
    return np.mean(ratios)  














'''
    SJS.
    The below class implements a module to compute dN/dS (or dN, dS individually) from mutation-selection model parameters, as described in Spielman and Wilke 2015 (doi: 10.1093/molbev/msv003).
'''


class dNdS_from_MutSel():
    '''
        Class to compute dN, dS, and/or dN/dS from mutation-selection model parameters, as described in Spielman and Wilke 2015.

        Positional arguments:
            codon_frequencies = either a list, numpy array, or dictionary of equilibrium codon frequencies. If list or numpy array, frequencies should be in alphabetical order, excluding stop codons. Relatively little sanity checking is done here, so provide something reasonable..
            mutation_dictionary = (optional) dictionary of mutation rates. This argument is analogous to the "mu" dictionary provided to pyvolve when simulating with custom mutation rates. If this argument is not provided, assumes equal mutation rates.
    '''
    def __init__(self, codon_frequencies, mutation_dictionary = None):
        
        self.ZERO = 1e-10
        self.codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}

        # Frequency setup
        self.codon_freqs = codon_frequencies
        if type(codon_frequencies) is dict:
            self.codon_freqs_dict = codon_frequencies
        else:
            self.codon_freqs_dict = {}
            for c in range(len(self.codons)):
                self.codon_freqs_dict[self.codons[c]] = self.codon_freqs[c]
        assert(abs(1. - np.sum(self.codon_freqs_dict.values())) < self.ZERO), "\n\nProvided codon frequencies must sum to 1."
        
        
        # Mutation setup        
        if mutation_dictionary is None:
            self.mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
        else:
            self.mu_dict = deepcopy(mutation_dictionary)
            for key in mutation_dictionary: 
                if key[1] + key[0] not in self.mu_dict:
                    self.mu_dict[key[1] + key[0]] = mutation_dictionary[key]



    def compute_dn(self):
        ''' 
            Compute dN from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        return self.compute_quantity("nonsyn")




    def compute_ds(self):
        ''' 
            Compute dS from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        return self.compute_quantity("syn")




    def compute_dnds(self):
        ''' 
            Compute dN/dS from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        dn = self.compute_dn()
        ds = self.compute_ds()
        
        if dn <= self.ZERO:
            return 0.
        elif ds <= self.ZERO:
            return np.inf
        else:
            return dn/ds


    def compute_quantity(self, type):
        '''
            Compute either dN or dS (type is 'nonsyn' or 'syn', respectively).
        '''
        numer = 0.
        denom = 0.

        for codon in self.codon_freqs_dict:
            rate = 0.
            sites = 0.
            rate, sites = self.calc_paths(codon, type)
            numer += rate
            denom += sites
        
        assert( denom > self.ZERO ), "\n\nProvided frequencies indicate no evolution is 'possible'."
        return numer/denom




     
    def calc_paths(self, source, type):
        ''' 
            Compute a term in the quantity numerator for a given source codon.
        '''
        rate = 0.
        sites = 0.
        source_freq = self.codon_freqs_dict[source]
        for target in self.codons:
            diff = self.get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
            if (type == 'nonsyn' and self.codon_dict[source] != self.codon_dict[target]) or (type == 'syn' and self.codon_dict[source] == self.codon_dict[target]):
                if len(diff) == 2:
                    rate  += self.calc_subst_prob( source_freq, self.codon_freqs_dict[target], self.mu_dict[diff], self.mu_dict[diff[1]+diff[0]] )
                    sites += self.mu_dict[diff]
        rate  *= source_freq
        sites *= source_freq
        return rate, sites



    def get_nuc_diff(self, source, target):
        '''
            Determine nucleotide difference between two codons.
        '''
        diff = ''
        for i in range(3):
            if source[i] != target[i]: 
                diff += source[i]+target[i]
        return diff

 
    
    def calc_subst_prob(self, pi, pj, mu_ij, mu_ji):
        '''
            Compute the substitution probability between two codons, as in Halpern and Bruno 1998.
        '''
        if abs(pi) <= self.ZERO or abs(pj) <= self.ZERO:
            fixation_rate = 0.
        else:
            p_mu = (mu_ji*pj)/(mu_ij*pi)
        
            # If p_mu == 1, L'Hopitals gives fixation rate of 1 (substitution probability is the forward mutation rate) 
            if abs(1. - p_mu) <= self.ZERO:
                fixation_rate = 1. 
            else:
                fixation_rate =  (np.log(p_mu))/(1. - (1./p_mu))
        return fixation_rate * mu_ij            





