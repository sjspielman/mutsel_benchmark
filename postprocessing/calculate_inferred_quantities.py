"""
    SJS
    This script computes and saves various quantities from inferences. These include the following:
        1. JSD between true and inferred AA frequency distributions
        2. Sum of absolute differences between true and inferred AA frequency distributions
        3. dN/dS
        4. Entropy of inferred AA frequency distributions
        5. Scaled selection coefficient distributions
"""
 
import sys
import numpy as np
sys.path.append("../")
from universal_functions import *

  
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



def dnds_from_params(site_fitness, mu_dict):
    """
        Calculate dN/dS from inferred parameters.
    """
    
    # Build the MutSel matrix (assumes equal codon frequencies per amino acid) with these parameters and extract equilibrium codon frequencies
    # Note that if mutation is symmetric, this gives the same results as Boltzmann distribution does anyways
    eqfreqs = codon_freqs_from_fitness_eigenvector(site_fitness, mu_dict)
    
    # Derive dN/dS and return
    c = dNdS_from_MutSel( dict(zip(g.codons, eqfreqs)), mu_dict)
    return c.compute_dnds()



def calculate_sum_absolute_diff(a, b):
    '''
        Compute the sum of absolute differences between amino acid frequency distributions a and b.
    '''
    return np.sum( [abs(a[i]-b[i]) for i in range(len(a))])


  
def calculate_kl(a,b):
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
    term1 = 0.5*calculate_kl(p, m)
    term2 = 0.5*calculate_kl(q, m)

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
        
        # Renormalize (as needed for entropy and KL calculates, what with division) to replace 0's with 1e-10
        true_aa_freqs[true_aa_freqs <= ZERO] = ZERO
        inferred_aa_freqs[inferred_aa_freqs <= ZERO] = ZERO
        true_aa_freqs /= np.sum(true_aa_freqs)
        inferred_aa_freqs /= np.sum(inferred_aa_freqs)
        
        
        jsd.append( calculate_jsd(true_aa_freqs, inferred_aa_freqs) )     
        diffsum.append( calculate_sum_absolute_diff(true_aa_freqs, inferred_aa_freqs) )  
        dnds.append( dnds_from_params(fitness[i], mu_dict) )
        entropy.append( calculate_entropy(inferred_aa_freqs) )
           
    assert(len(jsd)     == len(true)), "\nJSD not fully calculated for a dataset."
    assert(len(diffsum) == len(true)), "\nSum of absolute differences between frequencies not fully calculated for a dataset."
    assert(len(dnds)    == len(true)), "\ndN/dS not fully calculated for a dataset."
    assert(len(entropy) == len(true)), "\nEntropy not fully calculated for a dataset."

    with open(outfile, "w") as outf:
        outf.write("site,jsd,abs_sum_differences,dnds,entropy\n")
        outf.write( "\n".join([str(i+1)+","+str(jsd[i])+","+str(diffsum[i])+","+str(dnds[i])+","+str(entropy[i]) for i in range(len(jsd))]) )



def main():
    
    method_suffixes = {"swmutsel":"_MLE.txt", "phylobayes":".aap"}
    
    true_directory = "../simulation/true_simulation_parameters/"
    output_directory  = "dataframes/"
    
    # Grab all the true codon frequencies
    true_frequencies = {} 
    dms_datasets = ["HA", "NP", "LAC", "GAL4"]
    all_datasets = ["HA", "NP", "LAC", "GAL4", "1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A"]
    for data in all_datasets:
        t = np.loadtxt(true_directory + data + "_true_codon_frequencies.txt")
        true_frequencies[data] = t            

    for method in method_suffixes:
        input_directory = "../results/" + method + "/"
        suffix = method_suffixes[method] 
        filenames = os.listdir(input_directory)         
        for file in filenames:
            prefix = None
            if file.endswith(suffix):
                prefix = file.split(suffix)[0]
                is_dms = bool(len([x for x in range(len(dms_datasets)) if prefix.startswith(dms_datasets[x])]))
                if is_dms:
                    simprefix = prefix.split("_")[0]
                else:
                    simprefix = "_".join( prefix.split("_")[:2] )  
                print prefix, simprefix
                outfile_values = output_directory + prefix + "_statistics.csv"
                outfile_coeffs = output_directory + prefix + "_selcoeffs.csv"

                fitness, mu_dict = extract_parameters(input_directory, prefix)                         
                true_codon_frequencies = true_frequencies[simprefix]

                if not os.path.exists(outfile_values):
                    calculate_save_most_quantities(fitness, mu_dict, true_codon_frequencies, outfile_values)
                    
                #if not os.path.exists(outfile_coeffs):
                #    calculate_save_coeffs(fitness, outfile_coeffs)        

main()

    
    
