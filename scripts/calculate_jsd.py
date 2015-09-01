# SJS
# Script to compute the Jensen-Shannon Distance (JSD) between true and inferred fitness distributions for all simulated datasets.
# Before calculations, fitness values are re-scaled to be on the same linear scale as swMutSel, as this inference method fixes one of the fitnesses to 0. 

#usage = "Usage: python calculate_jsd.py <dataset> , where dataset is the name of a simulated dataset."
#assert(len(sys.argv) == 2), usage


from numpy import *
import sys
sys.path.append("universal_functions/")
from parsing_functions import *
from dnds_mutsel_functions import *
from pyvolve import Genetics
g = Genetics()
ZERO = 1e-10 # NOTE: everything with a frequency <= this value will be replaced with ZERO, for accurate comparisons



        
def fitness_to_bolzfreq(codon_fitness):
    
    # Convert codon fitness to codon frequencies with Boltzmann.
    codon_freqs = np.zeros(61)
    count = 0
    for codon in g.codons:
        codon_freqs[count] = np.exp( codon_fitness[codon] )
        count += 1
    codon_freqs /= np.sum(codon_freqs)                   
    assert( abs(1. - np.sum(codon_freqs)) <= ZERO), "codon_freq doesn't sum to 1 in codon_fitness_to_freqs"
    return codon_freqs 
    


def obtain_swmutsel_freqs(fitfile, hyphyfile):

    swmutsel_aafreqs = []
    sitewise_fitness = np.loadtxt(fitfile)
    pi, kappa, mu_dict = extract_optimized_params(hyphyfile, return_tree = False)       
    for site_fitness in sitewise_fitness:
          eqfreqs = codon_freqs_from_fitness(site_fitness, mu_dict)
          aafreqs = np.array(codon_to_aa_freqs(eqfreqs))
          aafreqs[aafreqs <= ZERO] = ZERO
          swmutsel_aafreqs.append(aafreqs)
    return swmutsel_aafreqs

# 
# def extract_true_fitness(file):
#     '''
#         Convert true codon frequencies to amino-acid fitness values.
#     '''
#     full_fitness = []
#     freqs = np.loadtxt(file)
#     for row in freqs:
#         row[row <= ZERO] = ZERO # replace 0 with ZERO (1e-10), since log on next line
#         cfit = np.log(row)
#     
#         aafit = []
#         for syn in g.genetic_code:
#             synindex = g.codons.index(syn[0])
#             aafit.append( cfit[synindex] )  
#         full_fitness.append(aafit) 
#         
#     return full_fitness 
#     
    
    
def calc_entropy(f):
    ''' 
        f is list of aa frequencies. Return Shannon entropy.
    '''
    return -1. * np.sum ( f[f > ZERO] * np.log(f[f > ZERO]) )   
    
    
    
def compute_sitewise_jsd(p1, p2):
    '''
        p1,p2 are each 2D arrays of site-specific AA frequency distributions.
        For each row in p1,p2, calculate the Jensen-Shannon Distance. Return list of sitewise distances.
    '''
    assert( len(p1) == len(p2) ), "Cannot compute JSD between sites in these two lists because the lists are different lengths."
    jsds = []
    for i in range(len(p1)):
        a = p1[i]
        b = p2[i]
        
        jsd = 0.
        for x in range(len(a)):
            stuff = (a[x] + b[x])/2.
            jsd += 0.5 * (a[x]*np.log(a[x]/stuff))
            jsd += 0.5 * (b[x]*np.log(b[x]/stuff))
        jsds.append( np.sqrt(jsd) )            
    return jsds
        
    
    
 
    

datasets = ["1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1RII_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated"]
outfile = "postprocessing/jsd_symmu.csv"
with open(outfile, "w") as f:
    f.write("dataset,site,jsd_true_sw\n")   #,jsd_true_pb\n")


for dataset in datasets:
    dataset_clean = dataset.replace("_simulated", "")

    # File names to read in
    truefile  = "simulation/flib/" + dataset_clean + "_codon_freq_lib.txt"
    swfitfile = "../results/simulation/" + dataset + "_nopenal_fitness.txt"
    swhyfile  = "../results/simulation/" + dataset + "_swmutsel_hyout.txt"
    #pbfitfile   = "../results/simulation/" + dataset + "_phylobayes.aap"

    # Obtain true amino-acid frequencies
    truefreqs = []
    rawtrue = np.loadtxt(truefile)
    for row in rawtrue:
        rawfreqs = np.array( codon_to_aa_freqs(row) )
        rawfreqs[rawfreqs <= ZERO] = ZERO
        truefreqs.append( rawfreqs )
    
    # Obtain swmutsel amino-acid frequencies
    #swfreqs = obtain_swmutsel_freqs(swfitfile, swhyfile)
    swfreqs = []
    rawsw = np.loadtxt(swfitfile)
    for row in rawsw:
        cfit = dict(zip(g.codons, aa_to_codon_fitness(row)))
        rawfreqs = fitness_to_bolzfreq(cfit)
        rawfreqs[rawfreqs <= ZERO] = ZERO
        swfreqs.append( rawfreqs )



    # Obtain pbmusel amino-acid frequencies
    ## code ##
    ## code ##
    
    
    # Calculate JSD between true/sw and true/pb
    true_sw_jsd = compute_sitewise_jsd(truefreqs, swfreqs)
    #true_pb_jsd = compute_sitewise_jsd(true_fitness, pb_fitness)
    

    print "Saving", dataset
    with open(outfile, "a") as f:
        for i in range(len(true_sw_jsd)):
            line = dataset + "," + str(i) + "," + str(true_sw_jsd[i]) + "\n" # "," + str(true_pb_jsd[i]) + "\n"
            f.write(line)
    
    
    
    
    
    
    
    
    
    
    
    
# 
# 
#     # Raw fitness values
#     true_fitness = extract_true_fitness(truefile)
#     sw_fitness   = np.loadtxt(swfile)
#     #pb_fitness   = np.loadtxt(pbfile)
#     
#     # Rescale true_fitness and pb_fitness to the sw_fitness linear scale.
#     for i in range(len(sw_fitness)):
#         j = np.where(sw_fitness[i] == 0.)[0]
#         true_fitness[i][j] = 0.
#         #pb_fitness[i][j]   = 0.
#     
# 
#     # Calculate JSD between true/sw and true/pb
#     true_sw_jsd = compute_sitewise_jsd(true_fitness, sw_fitness)
#     #true_pb_jsd = compute_sitewise_jsd(true_fitness, pb_fitness)
# 
# 
#     print "Saving", dataset
#     with open(outfile, "a") as f:
#         for i in range(len(true_sw_jsd)):
#             line = dataset + "," + str(i) + "," + str(true_sw_jsd[i]) + "\n" # "," + str(true_pb_jsd[i]) + "\n"
#             f.write(line)
#     assert 1==5
#             
# 
# 




