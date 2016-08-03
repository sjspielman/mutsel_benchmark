"""
     SJS
     This file contains functions which are generally used throughout this repository.
"""


import numpy as np
from pyvolve import *
from copy import deepcopy
g = Genetics()
ZERO=1e-10



def calculate_save_coeffs(fitness, outfile):
    '''
        Compute and save distribution of selection coefficients from a given set of fitness values.
    '''
    binned = []
    for site in fitness:
        for i in range(len(site)):
            f_i = site[i]
            for j in range(len(site)):
                diff = [x for x in range(3) if g.codons[i][x] != g.codons[j][x]]
                if len(diff) == 1 and i != j:
                    s = f_i - site[j]
                    if s <= -10:
                        s = -10.
                    if s >= 10.:
                        s = 10.
                    binned.append(s)
    with open(outfile, "w") as outf:
        outf.write("binnedcoeff\n")
        outf.write("\n".join([str(binned[x]) for x in range(len(binned))]) )




def calculate_entropy(f):
    '''
        Compute entropy of frequency distribution f.
    '''
    return -1. * np.sum ( f * np.log(f) )




def save_simulation_info(name, frequencies, fitnesses, omegas, entropies):
    '''
        Save simulation parameters to files.
    '''
    #freqfile = name + "_true_codon_frequencies.txt"
    #fitfile  = name + "_true_aa_fitness.txt"
    selcfile = name + "_true_selcoeffs.csv"
    #dnds_entropy_file = name + "_true_dnds_entropy.csv"

    #np.savetxt(freqfile, frequencies)
    #np.savetxt(fitfile, fitnesses)
    calculate_save_coeffs(fitnesses, selcfile)

    #with open(dnds_entropy_file, "w") as f:
    #    f.write("site,dnds,entropy\n")
    #    for i in range(len(omegas)):
    #        f.write(str(i+1)+","+str(omegas[i])+","+str(entropies[i])+"\n")



def codon_freqs_from_fitness_boltzmann(codon_fitness):
    '''
        Convert codon fitness values to stationary codon frequencies using Sella and Hirsh 2005 (Boltzmann).
        ***Assumes symmetric mutation rates***
    '''
    if type(codon_fitness) is not dict:
        codon_fitness = dict(zip(g.codons, codon_fitness))
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
        Extract equilibrium frequencies from eigenvector of MutSel matrix, built from fitness values and mutation rates.
    '''
    params = {"fitness": fitness, "mu": mu}
    m = Model("mutsel", params)
    return m.extract_state_freqs()


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
        aa_freqs.append(total)
    aa_freqs = np.array(aa_freqs)
    assert((1. - np.sum(aa_freqs)) <= ZERO), "\nImproperly converted codon to amino acid frequencies."
    return aa_freqs




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





#
# ########## Parsing functions for swMutSel and pbMutSel inferences, specifically for mutation rates ###########
#
#
# def compute_asym(mudict):
#     '''
#         Compute average mu_xy / mu_yx
#     '''
#     x = 0
#     completed = []
#     ratios = np.zeros(6)
#     for source in ["A", "C", "G", "T"]:
#         for target in ["A", "C", "G", "T"]:
#             if source == target or source+target in completed or target+source in completed:
#                 continue
#             else:
#                 ratio = mudict[source + target] / mudict[target + source]
#                 if ratio >= 1.:
#                     completed.append(source + target)
#                     ratios[x] = ratio
#                 else:
#                     completed.append(target + source)
#                     ratios[x] = 1./ratio
#                 x+=1
#     return np.mean(ratios)
