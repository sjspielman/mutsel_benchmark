# SJS
# Given a coding-sequence alignment and phylogeny, run the PhyloBayes mutsel model (Rodrigue 2010, 2014). 
# Once site-wise fitness values are inferred, compute site-wise dN/dS values.

# Usage: python run_phylobayes_mutsel.py <alignment_file> <tree_file> <cpu> <job_name> . 
## Alignment must be in phylip format, and tree must be in newick format. Both must be in working directory.
## Executables pb_mpi and readpb_mpi must be in the working directory.


import sys
import subprocess
from dendropy import Tree ## VERSION 4!!!!
 

def main():
    
    cpu = sys.argv[1]
    job_name = sys.argv[2]
    
    try:
        alnfile = sys.argv[3]
    except:
        assert(restart is True), "Specified alignment file does not exist. Path?"
    try:
        treefile = sys.argv[4]
    except:
        assert(restart is True), "Specified tree file does not exist. Path?"
    
    # Rewrite tree to create trifurcating root, as needed by phylobayes mpi
    tree = Tree.get_from_path(treefile, "newick", rooting = "force-unrooted")
    tree.resolve_polytomies() # in case of polytomies.
    tree.update_bipartitions() # this will create a trifurcating root on an unrooted tree
    tstring = str(tree).replace('[&U] ', '')
    with open('temp.tre', 'w') as tf:
        tf.write(tstring + ';\n')
        
    # Phylobayes is run to chain length 5500, sampling every 5 to yield 1100. Later, burnin of 100 is removed to get a final posterior n=1000 (same procedure as Rodrigue 2013 Genetics)
    pb_call = "mpirun -np " + str(cpu) + " ./pb_mpi -mutsel -cat -d " + alnfile + " -T temp.tre -x 5 1100 " + job_name
    
    run_pb_call = subprocess.call(pb_call, shell = True)
    assert( run_pb_call == 0 ), "pb_mpi didn't run!"
    
    # Parse output with readpb_mpi, using a burnin of 100 and saving everything else (posterior size = 1000)
    readpb_call = "mpirun -np " + str(cpu) + " ./readpb_mpi -x 100 1 -1 " + job_name + "\n"
    run_readpb_call = subprocess.call(readpb_call, shell = True)
    assert( run_readpb_call == 0 ), "readpb_mpi didn't run!"
    
    
main()
    


