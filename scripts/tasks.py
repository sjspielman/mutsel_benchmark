# Various tiny things we need to use Biopython for.

from Bio import AlignIO
import sys

file = sys.argv[1]
task = sys.argv[2]
# task == phy2fas: convert phylip to fasta
# task == alnlen: return alignment length from FASTA file.

if task == 'phy2fas':
    aln = AlignIO.read(file, 'phylip-relaxed')
    AlignIO.write(aln, file, 'fasta')

elif task == 'alnlen':
    aln = AlignIO.read(file, 'fasta')
    sys.stdout.write( str(len(aln[0])) )
