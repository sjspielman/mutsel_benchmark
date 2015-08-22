# SJS
# Given a coding-sequence alignment and phylogeny, run the swmutsel mutation-selection model (Tamuri 2012,2014) across a variety of penalty functions.
# Can either use HyPhy to optimize mutational parameters and branch lengths before sending to swmutsel, or just optimize branch lengths according to fixed mutation rates.
# Note that optimization is done first with HKY85 to get kappa and nucleotide frequencies, and then with a global MutSel neutral model (this is essentially MG94 with w=1), as this is how swmutsel interprets BL.
# Once site-wise fitness values are inferred, compute site-wise dN/dS values.

# Usage: python run_mutsel.py <alignment_file> <tree_file> <cpu> <dataset> <estimate_mutation>. 
## Alignment must be in either phylip/fasta format, and tree must be in newick format. Both must be in working directory.
## Additionally note that all executables called (swmutsel.jar and HYPHYMP) must be in the working directory.


import re
import os
import sys
import subprocess
import numpy as np
from Bio import AlignIO
from dnds_mutsel_functions import *

# Penalty functions (no penalty + 3 dirichlet + 3 multivariate normal) to run swmutsel with
penalty = {"nopenal" :"", 
           "d0.01"   :" -p dirichlet,0.01", 
           "d0.1"    :" -p dirichlet,0.1", 
           "d1.0"    :" -p dirichlet,1.0", 
           "mvn10"   :" -p mvn,10", 
           "mvn100"  :" -p mvn,100", 
           "mvn1000" :" -p mvn,1000"
          }

    
def optimize_tree_only(cpu, pi, kappa):
    '''
        Call hyphy batchfile "optimize_fmutsel_neutral_tree_only.bf" to obtain optimized branch lengths for *fixed nucleotide frequencies and kappa*
        Parse resulting output, and return all parameters (a new newick tree string, freqs, kappa). Note that freqs are returned as an ordered (T,C,A,G) csv string and kappa is also a string.
        NOTE: Assumes that the executable HYPHYMP is in the current directory.
    '''
    
    # Input the desired parameters to fix
    with open("optimize_fmutsel_neutral_tree_only.bf", "r") as f:
        raw = f.read()
    raw = raw.replace("insert_kappa", str(kappa))
    keypi = dict(zip(["insert_fA", "insert_fC", "insert_fG", "insert_fT"], pi))
    for key in keypi:
        raw = raw.replace(key, str(keypi[key]))
    with open("optimize_fmutsel_neutral_tree_only.bf", "w") as f:
        f.write(raw)
        
    # Optimize branch lengths
    run_hyphy = subprocess.call("./HYPHYMP CPU=" + str(cpu) + " optimize_fmutsel_neutral_tree_only.bf > hyout.txt", shell=True)
    assert(run_hyphy == 0), "optimize_fmutsel_neutral_tree_only.bf did not run!"

    with open("hyout.txt", "r") as hyout:
        hyout_string = str( hyout.read() )
        
    # Tree w/ optimized branch lengths #
    find_tree = re.search("Tree codon_tree=(.+)$", hyout_string)
    tree = find_tree.group(1)
    tree = re.sub(r"Node\d+", r"", tree) # remove node labels
    #tree = re.sub(r"([XN]P)(\d+)(\d)([XN]M)(\d+)(\d)", r"\1_\2.\3_\4_\5.\6", tree) # add back underscores, periods to taxon labels. This is necessary for the amine dataset and shouldn't match anything else.
    
    return tree
    
    
    

def optimize_fmutsel_neutral(cpu):
    '''
        Call hyphy batchfile "optimize_fmutsel_neutral.bf" to obtain optimized branch lengths, nucleotide frequencies, and kappa.
        Parse resulting output, and return all parameters (a new newick tree string, freqs, kappa). Note that freqs are returned as an ordered (T,C,A,G) csv string and kappa is also a string.
        NOTE: Assumes that the executable HYPHYMP is in the current directory.
    '''
    run_hyphy = subprocess.call("./HYPHYMP CPU=" + str(cpu) + " optimize_fmutsel_neutral.bf > hyout.txt", shell=True)
    assert(run_hyphy == 0), "optimize_fmutsel_neutral.bf did not run!"
    
    with open("hyout.txt", "r") as hyout:
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
    
    # Tree w/ optimized branch lengths #
    find_tree = re.search("Tree codon_tree=(.+)$", hyout_string)
    tree = find_tree.group(1)
    tree = re.sub(r"Node\d+", r"", tree) # remove node labels
    #tree = re.sub(r"([XN]P)(\d+)(\d)([XN]M)(\d+)(\d)", r"\1_\2.\3_\4_\5.\6", tree) # add back underscores, periods to taxon labels. This is necessary for the amine dataset and shouldn't match anything else.

    return nuc_freqs_ordered, kappa, tree, mu_dict



def build_mu_dict(pi, kappa):
    ''' Construct mutation rates dictionary.'''
    a = pi[0]
    c = pi[1]
    g = pi[2]
    t = pi[3]
    mu = {'AG':kappa*g, 'TC':kappa*c, 'GA':kappa*a, 'CT':kappa*t, 'AC':c, 'TG':g, 'CA':a, 'GT':t, 'AT':t, 'TA':a, 'GC':c, 'CG':g}  
    return mu
    

def run_swmutsel(pi, kappa, alnfile, treefile, cpu, dataset, penalname, penalarg, label):
    '''
        Call swmutsel and parse output. Return site-wise amino acid fitnesses as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ]
        By default, uses the Dirichlet penalty function with shape 0.1.
        NOTE: Assumes that swmutsel executable is in the current directory and named "swmutsel.jar"   
    '''
    job_name = dataset + '_' + penalname + label
        
    runstring = "java -jar swmutsel.jar -n " + job_name + " -s " + alnfile + " -t " + treefile + " -fix mutation,branches -k " + kappa + " -pi " + pi + " -gc standard "  + " -threads " + str(cpu) + penalarg
    run_swmutsel = subprocess.call(runstring, shell=True)
    assert(run_swmutsel == 0), "swmutsel did not run!"
    
    fitness = []
    new_order = [0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17] # need to reorder fitnesses from tamuri's mapping to mine.
    with open(job_name + '_MLE.txt', 'r') as f:
        for line in f:
            newline = line.split(',')[1:] # first col in csv is site index, so ignore it.
            if len(newline) == 20:
                fitness.append( [float(y) for (x,y) in sorted(zip(new_order,newline))] )
    
    # Save a file with fitness values
    np.savetxt(job_name + '_fitness.txt', fitness, delimiter = '\t')        
    
    return fitness    









def main():
    usage = "\nUsage: python run_swmutsel.py <dataset> <cpu> <estimate_mutation>. Note that all files and all executables ('swmutsel'+'HYPHYMP') must be in the working directory!"
    assert( len(sys.argv) == 4 ), usage
    
    dataset = sys.argv[1]
    cpu = sys.argv[2]
    estimate_mutation = bool(int(sys.argv[3]))

    alnfile_fasta = dataset + ".fasta"
    alnfile_phy = dataset + ".phy"
    treefile = dataset + ".tre"
    #assert( os.path.exists(alnfile_fasta) ), "There is no fasta alignment file."
    #assert( os.path.exists(alnfile_phy) ), "There is no phylip alignment file."
    #assert( os.path.exists(treefile) ), "There is tree file."

    if estimate_mutation is True:
        label = "_estmu"
    else:
        label = "_fixedmu"
    opt_treefile = dataset + label + "_optimized.tre"


    
    
    # Prep hyphy input file and call hyphy to optimize mutational parameters, create mu_dict to use later for dnds derivation, and make a treefile with the optimized tree
    subprocess.call("cat " + alnfile_fasta + " " + treefile + " > hyin.txt", shell=True)

    
    if estimate_mutation is True:
        pi, kappa, treestring, mu_dict = optimize_fmutsel_neutral(cpu)    
    
    else:
        pi = ['0.25', '0.25', '0.25', '0.25']
        kappa = 4.0; mu = 1.
        mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}
        treestring = optimize_tree_only(cpu, pi, kappa)
   
    pi2 = [str(i) for i in pi]
    pi_string = ",".join(pi2)             
    with open(opt_treefile, "w") as f:
        f.write(treestring)



    # Call swmutsel to obtain site-wise fitness values across a variety of penalizations (no penalty and those tested in Tamuri 2014).
    for penal in penalty:
        
        sitewise_fitness = run_swmutsel( pi_string, str(kappa), alnfile_phy, opt_treefile, cpu, dataset, penal, penalty[penal], label )
    
        # Compute site-wise dN/dS from amino-acid fitness values
        dnds = []
        for site_fitness in sitewise_fitness:
            
            # Build the MutSel model with these parameters and extract equilibrium codon frequencies
            codon_fitness = aa_to_codon_fitness(site_fitness)
            matrix = build_mutsel_matrix(mu_dict, codon_fitness)
            eqfreqs = get_eq_from_eig(matrix)
            cf = dict(zip(codons, eqfreqs))
            
            # Derive dN/dS and append to list
            dnds.append( derive_dnds(cf, mu_dict)  )    

        # Save
        with open(dataset + "_" + penal + label + "_dnds.txt", "w") as outf:
            outf.write( "\n".join([str(i) for i in dnds]) )
  
    
main()
    


