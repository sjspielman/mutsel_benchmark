# SJS
# Script to compute statistics about empirical datasets. 

from dendropy import Tree
from dendropy.calculate import treemeasure # use this line for PD
from Bio import AlignIO
import numpy as np
from universal_functions import *

datadir = "../data/empirical/"
rawdir  = "../results/raw_results/empirical/"
outfile = "../results/empirical_data_statistics.csv"
with open(outfile, "w") as f:
    f.write("dataset,ntaxa,ncol,treelen,meanpd,asym_nopenal,gc_nopenal\n")
    
    
datasets = ["PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715", "pb2"]
for data in datasets:
    
    print data
    
    ##### Alignment information #####
    aln = AlignIO.read(datadir + data + ".fasta", "fasta")
    ntaxa = str(len(aln))
    ncol  = str(len(aln[0]))
    
    ##### Tree information #####
    t = Tree.get_from_path(datadir + data + ".tre", 'newick')
    tree_length = str(t.length())
    
    pd = []
    dist = treemeasure.PatristicDistanceMatrix(tree=t)
    for i, t1 in enumerate(t.taxon_namespace):
        for t2 in t.taxon_namespace[i:]:
            d = dist(t1,t2)
            pd.append( float(d) )
    meanpd = str( np.mean(pd) )
       
        
    
    ##### Mutational asymmetry #####
    
    # swMutSel, no penalty
    mu_dict = parse_swMutSel_mutation(rawdir + "swmutsel/" + data + "_nopenal_MLE.txt")   
    asym = str(compute_asym(mu_dict))
    
# 
#     # swMutSel, d1.0
#     mu_dict = parse_swMutSel_mutation(rawdir + "swmutsel/" + data + "_d1.0_MLE.txt")   
#     d_asym = str(compute_asym(mu_dict))
#         
#     # pbMutSel
#     try:
#         mu_dict = parse_pbMutSel_mutation(rawdir + "phylobayes/" + data + "_phylobayes.trace")
#         pb_asym = str(compute_asym(mu_dict))
#     except:
#         pb_asym = "NA"
#     
    
    ##### Base composition (as GC) #####
    
    # swMutSel, nopenalty
    with open(rawdir + "swmutsel/" + data + "_nopenal_MLE.txt", "r") as f:
        lines = f.readlines()
    rawpis   = lines[13].strip().split(",")
    pis = []
    for pi in rawpis:
        pis.append(float(pi)) 
    gc = str(pis[1] + pis[3])
    
    

    ##### Save #####
    save_line = data + "," + ntaxa + "," + ncol + "," + tree_length + "," + meanpd + "," + asym + "," + gc + "\n"
    with open(outfile, "a") as f:
        f.write(save_line)

    
    
    
    
        