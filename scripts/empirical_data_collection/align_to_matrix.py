# SJS 2/13/15.
# Script to clean up alignments:
### 1. Remove columns with >=1% missing/ambiguous data
### 2. Remove rows with >=1% missing/ambiguous data
### 3. Remove duplicate sequences (using external call to RAxML in script dup_hack.sh to do it for me, hurray!!)
### 4. Save to ../../

from Bio import AlignIO
import numpy as np
import sys,os,subprocess
from string import maketrans
np.set_printoptions(threshold='nan')        

class AlignMatrix(object):
    
    
    def __init__(self, **kwargs):
        
        alnfile  = kwargs.get("alnfile", "")     
        alnfmt   = kwargs.get("format", "fasta")
        rawcode  = kwargs.get("type", "").lower() # dna, protein, codon

        
        self.code   = self._set_code(rawcode)
        self.step   = len(self.code.keys()[0])
        self.rawaln = AlignIO.read(alnfile, alnfmt)
        self.numseq = len(self.rawaln)
        self.alnlen = len(self.rawaln[0]) / self.step
        self.id_dict = {}
        self.aln_dict = {}


    def alnmat_numseq(self):
        return self.alnmat.shape[0]

    def alnmat_alnlen(self):
        return self.alnmat.shape[1]
    


    def _set_code(self, type):
        '''
            Gaps are given 0!
        '''
        dna     = ["-", "A", "C", "G", "T"]
        protein = ["-", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        codon   = ["---", "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        l = eval(type)
        
        return dict(zip( l, np.arange(len(l)) )) 
        

    def seq_to_int(self, row, id_count):
        '''
            Replace sequence with integer equivalent. All gaps are 0, and all missing are negatives, and self.code is built up accordingly.
            Again note, the first entry in a given row corresponds to the sequence ID.
        '''
        new_row = np.zeros(self.alnlen + 1)
        new_row[0] = id_count
        for i in range(0, len(row), self.step):
            pos = row[i:i+self.step]
            try:
                value = self.code[pos]
            except:
                value = min(self.code.values()) - 1
                self.code[pos] = value
            new_row[(i/self.step) + 1] = value
        return new_row


    def int_to_seq(self, mapdict, row):
        '''
            Replace integer with actual sequence.
        '''
        new_row = ""
        for pos in row:
            new_row += mapdict[pos]
        return new_row
        
        

    def convert_aln_to_matrix(self):
        '''
            Convert biopython AlignIO object into numpy integer matrix, where rows are sequences and columns are, indeed, columns.
            Note that the first column contains an integer mapping to the sequence id.
        '''
        
 
        self.alnmat = np.zeros( [ self.numseq, self.alnlen + 1 ], dtype = 'int16')
        i = 0
        for record in self.rawaln:
            id = str(record.id)
            self.id_dict[i] = id
            seq = str(record.seq)
            self.alnmat[i] = self.seq_to_int( seq, i )
            #print self.alnmat[i][0], i
            i += 1



    def convert_matrix_to_aln(self):
        '''
            From an alignment matrix, create an alignment dictionary with original ids as keys and updated sequences as values.
        '''
        inv_code = {v: k for k, v in self.code.items()}
        for row in self.alnmat:
            id = self.id_dict[ row[0] ]
            self.aln_dict[ id ] = self.int_to_seq( inv_code, row[1:] )
        

        
      
      
      
      
      
        
        
        
        
        
        