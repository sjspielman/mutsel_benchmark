# SJS
# Extract codon frequencies, dnds from deep mutational scanning AA preferences combined with experimental mutation rates. 
# Note that these parameters, as they are experimental, are simply taken as they are, and NOT subjected to a strong/weak regime.


import sys
sys.path.append("../")
from universal_functions import *
import numpy as np
from scipy import linalg
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codon_dict  = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
codons      = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
ZERO        = 1e-10


def get_nuc_diff(source, target, grab_position = False):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff
        
        
        
def build_metropolis_matrix(amino_prop_dict, mu_dict):
    ''' metropolis only, as this matrix definition more suits experimental propensities according to Bloom 2014. '''
    matrix = np.zeros([61,61])
    
    # off-diagonal entries
    for x in range(61):
        x_codon = codons[x]
        x_aa = codon_dict[x_codon]
        fx = amino_prop_dict[x_aa]
        for y in range(61):
            y_codon = codons[y]
            y_aa = codon_dict[y_codon]
            fy = amino_prop_dict[y_aa]
            diff = get_nuc_diff(x_codon, y_codon)
            if len(diff)==2:
                if x_aa == y_aa or fy >= fx:
                    matrix[x][y] =  mu_dict[diff]
                else:
                    matrix[x][y] = fy / fx  * mu_dict[diff]    
    # diagonal entries
    for i in range(61):
        matrix[i][i] = -1. * np.sum(matrix[i]) 
        assert( -ZERO < np.sum(matrix[i]) < ZERO ), "diagonal fail"
    return matrix


def get_eq_from_eig(m):   
    ''' get the equilibrium frequencies from the matrix. the eq freqs are the *left* eigenvector corresponding to eigenvalue of 0. 
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
    max_v = max_v.real # these are the stationary frequencies
    assert( abs(np.sum(max_v) - 1.) < ZERO), "Eigenvector of equilibrium frequencies doesn't sum to 1."

    # SOME SANITY CHECKS
    assert np.allclose(np.zeros(61), np.dot(max_v, m)) # should be true since eigenvalue of zero
    pi_inv = np.diag(1.0 / max_v)
    s = np.dot(m, pi_inv)
    assert np.allclose(m, np.dot(s, np.diag(max_v)), atol=ZERO, rtol=1e-5), "exchangeability and equilibrium does not recover matrix"
    # additional overkill check
    for i in range(61):
        pi_i = max_v[i]
        for j in range(61):
            pi_j = max_v[j]
            forward  = pi_i * m[i][j] 
            backward = pi_j * m[j][i]
            assert(abs(forward - backward) < 1e-6), "Detailed balance violated."  # note that we need to use high tolerance here because the propensities have very few digits so we encounter more FLOP problems.
    return max_v


mudict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
truedir = "true_simulation_parameters/"

for source in ["HA", "NP"]:
    
    infile               = truedir + source + "_preferences.txt"
    outfile_freqs        = truedir + source + "_true_codon_frequencies.txt"
    outfile_dnds_entropy = truedir + source + "_true_dnds_entropy.csv"
    
    raw_preferences = np.loadtxt(infile)
    nsites = len(raw_preferences)

    final_codon_freqs = np.zeros([nsites, 61])
    final_dnds = np.zeros(nsites)
    final_entropy = np.zeros(nsites)
    for i in range(nsites):
        print i
        amino_prefs_dict = dict(zip(amino_acids, raw_preferences[i])) 
        m = build_metropolis_matrix(amino_prefs_dict, mudict)
        cf = get_eq_from_eig(m) 
        assert( abs(np.sum(cf) - 1.) < ZERO ), "codon frequencies do not sum to 1" 

        c = dNdS_from_MutSel(dict(zip(codons, cf)), mudict)
        dnds = c.compute_dnds()

        final_codon_freqs[i] = cf
        final_dnds[i] = dnds
        final_entropy[i] = calc_entropy( codon_freqs_to_aa_freqs(cf)) 

    np.savetxt(outfile_freqs, final_codon_freqs)
    with open(outfile_dnds_entropy, "w") as f:
        f.write("site,dnds,entropy\n")
        for x in range(len(final_dnds)):
            f.write(str(x+1) + "," + str(final_dnds[x]) + "," + str(final_entropy[x]) + "\n")
