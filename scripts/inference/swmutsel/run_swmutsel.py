
import re
import sys
import subprocess
from numpy import savetxt
from universal_functions import *


# Penalty functions (no penalty + 3 dirichlet + 3 multivariate normal) to run swmutsel with
penalty = {"nopenal" :"", 
           "d0.01"   :" -p dirichlet,0.01", 
           "d0.1"    :" -p dirichlet,0.1", 
           "d1.0"    :" -p dirichlet,1.0", 
           "mvn10"   :" -p mvn,10", 
           "mvn100"  :" -p mvn,100", 
           "mvn1000" :" -p mvn,1000"
          }
    

def run_swmutsel(pi, kappa, alnfile, treefile, cpu, dataset, penalname, penalarg):
    '''
        Call swmutsel and parse output. Return site-wise amino acid fitnesses as list: [ [site1_fitnesses], [site2_fitnesses], [site3_fitnesses] ... [siten_fitnesses] ]
        By default, uses the Dirichlet penalty function with shape 0.1.
        NOTE: Assumes that swmutsel executable is in the current directory and named "swmutsel.jar"   
    '''
    job_name = dataset + '_' + penalname
        
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
    
    # Save a file with the correctly-ordered (well, my order) amino-acid fitness values
    savetxt(job_name + '_fitness.txt', fitness, delimiter = '\t')        









def main():
    usage = "\nUsage: python run_swmutsel.py <aln> <treefile> <cpu> <type>. Note that all files and all executables ('swmutsel'+'HYPHYMP') must be in the working directory!"
    assert( len(sys.argv) == 5 ), usage
    
    dataset = sys.argv[1]
    treefile = sys.argv[2]
    cpu = sys.argv[3]
    type = sys.argv[4]
    
    if type == "emp":
        run_hyphy = True
    elif type == "sim":
        run_hyphy = False
    else:
        raise ValueError("type (4th argument) must be either 'emp' or 'sim'.")

    alnfile_fasta = dataset + ".fasta"
    alnfile_phy = dataset  + ".phy"
    
    if run_hyphy:
        # Prep hyphy input file and call hyphy to optimize mutational parameters, create mu_dict to use later for dnds derivation, and make a treefile with the optimized tree
        subprocess.call("cat " + alnfile_fasta + " " + treefile + " > hyin.txt", shell=True)
        run_hyphy = subprocess.call("./HYPHYMP CPU=" + str(cpu) + " optimize_fmutsel_neutral > hyout.txt", shell=True)
        assert(run_hyphy == 0), "Hyphy optimization did not run."
        
        pi, kappa, treestring, mu_dict = extract_optimized_params("hyout.txt")  
        pi2 = [str(i) for i in pi]
        pi_string = ",".join(pi2)             

        opt_treefile = dataset + "_swmutsel_optimized.tre"
        with open(opt_treefile, "w") as f:
            f.write(treestring)

    else:
        opt_treefile = treefile
        pi_string = "0.25,0.25,0.25,0.25"
        kappa = "1.0"
         

    # Call swmutsel to obtain site-wise fitness values across a variety of penalizations (no penalty and those tested in Tamuri 2014).
    for penal in penalty:
        run_swmutsel( pi_string, str(kappa), alnfile_phy, opt_treefile, cpu, dataset, penal, penalty[penal])
    
main()
    


