# Quick script to convert alignment from fasta to phylip

from Bio import AlignIO
import sys

infile = sys.argv[1]
outfile = infile.split(".fasta")[0] + ".phy"
AlignIO.convert(infile, "fasta", outfile, "phylip-relaxed")