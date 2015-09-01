# SJS
# Given a coding-sequence alignment and phylogeny, run the PhyloBayes mutsel model (Rodrigue 2010, 2014). 
# Once site-wise fitness values are inferred, compute site-wise dN/dS values.

# Usage: python run_phylobayes_mutsel.py <alignment_file> <tree_file> <cpu> <job_name> . 
## Alignment must be in phylip format, and tree must be in newick format. Both must be in working directory.
## Executables pb_mpi and readpb_mpi must be in the working directory.


import sys
from dendropy import Tree
 

def main():
    
    restart = bool(sys.argv[1])
    cpu = sys.argv[2]
    job_name = sys.argv[3]
    
    try:
        alnfile = sys.argv[4]
    except:
        assert(restart is True), "Specified alignment file does not exist. Path?"
    try:
        treefile = sys.argv[5]
    except:
        assert(restart is True), "Specified tree file does not exist. Path?"
    
    
    # Run phylobayes, either from scratch or restarting a previous chain
    # Phylobayes is run to chain length 5500, sampling every 5 to yield 1100. Later, burnin of 100 is removed to get a final posterior n=1000
    if not restart:
        # Rewrite tree to create trifurcating root, as needed by phylobayes mpi
        tree = Tree.get_from_path(treefile, "newick", rooted = False)
        tree.resolve_polytomies() # in case of polytomies.
        tree.update_splits() # this will create a trifurcating root on an unrooted tree
        tstring = str(tree).replace('[&U] ', '')
        with open('temp.tre', 'w') as tf:
            tf.write(tstring + ';\n')
            
        pb_call = "mpirun -np " + str(cpu) + " ./pb_mpi -mutsel -cat -d " + alnfile + " -T " + treefile + " -x " + str(every) + " " + str(until) + " " + job_name
    
    else:
        pb_call = "mpirun -np 16 " + str(cpu) + " ./pb_mpi " + job_name    
    
    run_pb_call = subprocess.call(pb_call, shell = True)
    assert( run_pb_call == 0 ), "pb_mpi didn't run!"
    
    # Parse output with readpb_mpi, using a burnin of 100 and saving everything else (posterior size = 1000)
    readpb_call = "mpirun -np " + str(cpu) + " ./readpb_mpi -x 100 1 -1 " + job_name + "\n"
    run_readpb_call = subprocess.call(readpb_call, shell = True)
    assert( run_readpb_call == 0 ), "readpb_mpi didn't run!"
    
    
main()
    


