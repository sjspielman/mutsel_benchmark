
ZERO = 1e-10
import re
import numpy as np
from scipy import linalg
import sys
from random import uniform, shuffle
from pyvolve import *
g = Genetics()


def aa_to_codon_fitness(fitness):
    '''
        Convert list of amino acid *fitnesses* to list of codon *fitness*
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


def aa_to_codon_freqs(aa_freqs):
    
    codon_freqs = np.zeros(61)
    aa_dict = dict(zip(g.amino_acids, aa_freqs))
    for aa in aa_dict:
        syn_codons = g.genetic_code[ g.amino_acids.index(aa) ]
        cf = aa_dict[aa] / float(len(syn_codons))
        for syn in syn_codons:
            codon_freqs[ g.codons.index(syn) ] = cf
    
    assert(1. - np.sum(codon_freqs) <= ZERO), "bad codon freqs from aa freqs"
    return codon_freqs, dict(zip(g.codons, codon_freqs))


def get_eq_from_eig(m):   
    ''' get the equilibrium frequencies from the matrix. the eq freqs are the left eigenvector corresponding to eigenvalue of 0. 
        Code here is largely taken from Bloom. See here - https://github.com/jbloom/phyloExpCM/blob/master/src/submatrix.py, specifically in the fxn StationaryStates
    '''
    (w, v) = linalg.eig(m, left=True, right=False)
    max_i = 0
    max_w = w[max_i]
    for i in range(1, len(w)):
        if w[i] > max_w:
            max_w = w[i]
            max_i = i
    assert( abs(max_w) < ZERO ), "Maximum eigenvalue is not close to zero."
    max_v = v[:,max_i]
    max_v /= np.sum(max_v)
    eq_freqs = max_v.real # these are the stationary frequencies
    
    # SOME SANITY CHECKS
    assert np.allclose(np.zeros(61), np.dot(eq_freqs, m)) # should be true since eigenvalue of zero
    pi_inv = np.diag(1.0 / eq_freqs)
    s = np.dot(m, pi_inv)
    assert np.allclose(m, np.dot(s, np.diag(eq_freqs)), atol=ZERO, rtol=1e-5), "exchangeability and equilibrium does not recover matrix"
    
    # And for some impressive overkill, double check pi_i*q_ij = pi_j*q_ji
    for i in range(61):
        pi_i = eq_freqs[i]
        for j in range(61):
            pi_j = eq_freqs[j]
            forward  = pi_i * m[i][j] 
            backward = pi_j * m[j][i]
            assert(abs(forward - backward) < ZERO), "Detailed balance violated."    
    return eq_freqs




def build_mutsel_matrix(mu_dict, codon_fitness):
    ''' 
        Build the mutsel model matrix according to,  q_ij = mu_ij * Sij / (1-e^(-Sij)).
        The resulting *left* eigenvector of this matrix gives the state frequencies.
    '''
    mb = mutSel_Matrix(self.params, "neutral")
    self.matrix = mb()
    
    matrix = np.zeros([61,61])
    for s in range(61):
        for t in range(61):
        
            nucdiff = get_nuc_diff(codons[s], codons[t])
            if len(nucdiff) != 2:
                continue
                
            # Non-diagonal
            sij = codon_fitness[t] - codon_fitness[s]  
            if abs(sij) < ZERO:
                rate = mu_dict[nucdiff]
            else:
                rate = mu_dict[nucdiff] *  (sij)/(1 - np.exp(-1.*sij))
            matrix[s][t] = rate
                
        # Fill in the diagonal position so the row sums to 0, but ensure it doesn't became -0
        matrix[s][s]= -1. * np.sum( matrix[s] )
        if matrix[s][s] == -0.:
            matrix[s][s] = 0.
        assert ( abs(np.sum(matrix[s])) < ZERO ), "Row in instantaneous matrix does not sum to 0."
    return matrix

    




def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff







