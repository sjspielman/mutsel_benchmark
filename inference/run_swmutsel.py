
import re
import sys
import subprocess
from numpy import savetxt


# Penalty functions (no penalty + 3 dirichlet + 3 multivariate normal) to run swmutsel with
penalty = {"d0.01"   :" -p dirichlet,0.01", 
           "d0.1"    :" -p dirichlet,0.1", 
           "mvn10"   :" -p mvn,10", 
           "mvn100"  :" -p mvn,100",
           "nopenal" :""}

def run_swmutsel(alnfile, treefile, cpu, dataset, penalname, penalarg):
    '''
        Call swmutsel and parse output. Return site-wise amino acid fitnesses as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ]
        swMutSel is run under a provided penalty function. The pi vector and kappa are optimized, but the tree is fixed.
        NOTE: Assumes that swmutsel executable is in the current directory and named "swmutsel.jar"   
    '''
    job_name = dataset + '_' + penalname
        
    runstring = "java -jar swmutsel.jar -n " + job_name + " -s " + alnfile + " -t " + treefile + " -fix branches -gc standard "  + " -threads " + str(cpu) + penalarg
    run_swmutsel = subprocess.call(runstring, shell=True)
    assert(run_swmutsel == 0), "swmutsel did not run!"
    
    fitness = []
    new_order = [0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17] # need to reorder fitnesses from tamuri's mapping to mine.
    with open(job_name + '_MLE.txt', 'r') as f:
        for line in f:
            newline = line.split(',')[1:] # first col in csv is site index, so ignore it.
            if len(newline) == 20:
                fitness.append( [float(y) for (x,y) in sorted(zip(new_order,newline))] )
    
    # Save a file with the correctly-ordered (well, my order) amino-acid fitness values
    savetxt(job_name + '_fitness.txt', fitness, delimiter = '\t')        









def main():
    usage = "\nUsage: python run_swmutsel.py <aln> <treefile> <penalty> <cpu>. Note that all files and all executables must be in the working directory!"
    assert( len(sys.argv) == 5 ), usage
    
    dataset  = sys.argv[1]
    treefile = sys.argv[2]
    penal    = sys.argv[3]
    cpu      = sys.argv[4]
    alnfile  = dataset  + ".phy"

    run_swmutsel( alnfile, treefile, cpu, dataset, penal, penalty[penal])
    
main()
    


