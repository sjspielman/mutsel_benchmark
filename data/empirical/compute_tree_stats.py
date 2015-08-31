# SJS
# Script to compute statistics about empirical datasets. Results are printed in the form of a LaTex table, for convenience.


from dendropy import Tree
from dendropy.calculate import treemeasure # use this line for PD
from Bio import AlignIO

sep = " & "
datasets = ["amine", "pb2", "PF00126", "PF01336", "PF01926", "PF04055", "PF00593", "PF07715", "PF01266", "PF02518"]
for data in datasets:

    table_line = data + sep

    # Alignment information
    aln = AlignIO.read(data + ".fasta", "fasta")
    ntaxa = str(len(aln))
    ncol  = str(len(aln[0]))
    table_line += ntaxa + sep + ncol + sep
    
    # Tree information
    t = Tree.get_from_path(data + ".tre", 'newick')
    tree_length = str(round(t.length(), 2))
    pd = []
    dist = treemeasure.PatristicDistanceMatrix(tree=t)
    for i, t1 in enumerate(t.taxon_namespace):
        for t2 in t.taxon_namespace[i:]:
            d = dist(t1,t2)
            pd.append( float(d) )
    mean_pairwise = str(round(sum(pd) / float(len(pd)),2))

    
    table_line += tree_length + sep + mean_pairwise + "  \\\\"
    
    print table_line
        