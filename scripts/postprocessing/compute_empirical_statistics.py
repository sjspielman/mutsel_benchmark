# SJS
# Script to compute statistics about empirical datasets. Results are saved to a csv and also printed to stdout in the form of a LaTex table, for manuscript convenience.


from dendropy import Tree
from dendropy.calculate import treemeasure # use this line for PD
from Bio import AlignIO

datadir = "../../data/empirical/"
outfile = "empirical_data_statistics.csv"
with open(outfile, "w") as f:
    f.write("dataset,ntaxa,ncol,treelen,meanpairwise\n")
    
    
datasets = ["amine", "pb2", "PF00126", "PF00593", "PF01336", "PF01266", "PF01926", "PF02518", "PF04055", "PF07715"]
for data in datasets:
    
    # Alignment information
    aln = AlignIO.read(datadir + data + ".fasta", "fasta")
    ntaxa = str(len(aln))
    ncol  = str(len(aln[0]))
    
    # Tree information
    t = Tree.get_from_path(datadir + data + ".tre", 'newick')
    tree_length = str(round(t.length(), 2))
    pd = []
    dist = treemeasure.PatristicDistanceMatrix(tree=t)
    for i, t1 in enumerate(t.taxon_namespace):
        for t2 in t.taxon_namespace[i:]:
            d = dist(t1,t2)
            pd.append( float(d) )
    mean_pairwise = str(round(sum(pd) / float(len(pd)),2))


    # Save
    save_line = data + "," + ntaxa + "," + ncol + "," + tree_length + "," + mean_pairwise + "\n"
    with open(outfile, "a") as f:
        f.write(save_line)
    

    # Latex table
    latex_line = save_line.replace(",", " & ").strip()
    print latex_line  + "  \\\\"

    
    
    
        