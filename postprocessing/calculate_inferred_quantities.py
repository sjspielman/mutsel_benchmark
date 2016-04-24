'''
    SJS
    This script computes and saves various quantities from inferences. These include the following:
        1. JSD between true and inferred AA frequency distributions
        2. Sum of absolute differences between true and inferred AA frequency distributions
        3. dN/dS
        4. Entropy of inferred AA frequency distributions
        5. Scaled selection coefficient distributions
'''
 
import sys
import numpy as np

sys.path.append("../")
from universal_functions import *

  
  
def calc_sum_absolute_diff(a, b):
    '''
        Compute the sum of absolute differences between amino acid frequency distributions a and b.
    '''
    return np.sum( [abs(a[i]-b[i]) for i in range(len(a))])


  
def calc_kl(a,b):
    '''
        Compute KL distance from amino acid frequency distributions b to a.
    '''
    return np.sum( a * np.log(a/b) ) 

    

def calculate_jsd(p, q):
    '''
        Calculate JSD between frequency distributions.
        p and q are input distributions to compare, each of length 20.
    '''

    m = (p+q)/2.
    term1 = 0.5*calc_kl(p, m)
    term2 = 0.5*calc_kl(q, m)

    return np.sqrt( term1 + term2 )


def calculate_save_most_quantities(fitness, mu_dict, true, outfile):
    '''
        Compute and save summary quantities for inferences. Quantities include...
            1. JSD between true and inferred AA frequency distributions
            2. Sum of absolute differences between true and inferred AA frequency distributions
            3. dN/dS
            4. Entropy for AA distribution 
    '''
    
    jsd     = []
    diffsum = []
    dnds    = []
    entropy = []
    for i in range(len(fitness)):
        true_aa_freqs = np.array(codon_freqs_to_aa_freqs(true[i]))
        inferred_codon_freqs = np.array( codon_freqs_from_fitness_eigenvector(fitness[i], mu_dict) )
        inferred_aa_freqs = codon_freqs_to_aa_freqs(inferred_codon_freqs)
        
        jsd.append( calculate_jsd(true_aa_freqs, inferred_aa_freqs) )     
        diffsum.append( calculate_calc_sum_absolute_diff(true_aa_freqs, inferred_aa_freqs) )  
        dnds.append( dnds_from_params(fitness[i], mu_dict) )
        entropy.append( calc_entropy(inferred_aa_freqs) )
           
    assert(len(jsd) == len(true)), "\nJSD not fully calculated for a dataset."
    assert(len(diffsum) == len(true)), "\nSum of absolute differences between frequencies not fully calculated for a dataset."
    assert(len(dnds) == len(true)), "\ndN/dS not fully calculated for a dataset."
    assert(len(entropy) == len(true)), "\nEntropy not fully calculated for a dataset."

    with open(outfile, "w") as outf:
        outf.write("site,jsd,abs_sum_differences,dnds,entropy\n")
        outf.write( "\n".join([str(i+1)+","+str(jsd[i])+","+str(diffsum[i])+","+str(dnds[i])+","+str(entropy[i]) for i in range(len(jsd))]) )



def main():
    
    method_suffixes = {"swmutsel":"_MLE.txt", "phylobayes":".aap"}
    
    truedir = "../simulation/true_simulation_parameters/"
    rawdir  = "../results/"
    outdir  = rawdir + "dnds_coeffs_jsd/"
    
    for method in method_suffixes:
        indir = rawdir + method + "/"
        suffix = method_suffixes[method]         
        filenames = os.listdir(input_directory)
        for file in filenames:
    
            prefix = None
            if file.endswith(suffix):
                prefix = file.split(suffix)[0]
            
                outfile_values = output_directory + prefix + "_dnds_entropy_jsd_diffsum.csv"
                outfile_coeffs = output_directory + prefix + "_selcoeffs.csv"

                fitness, mu_dict = extract_parameters(input_directory, prefix)                         
                true_codon_frequencies = np.loadtxt(true_directory + simprefix + "_true_codon_frequencies.txt")

                calculate_save_most_quantities(fitness, mu_dict, true_codon_frequencies, outfile_values)
                calculate_save_coeffs(fitness, outfile_coeffs)        

main()

    
    
