# This script simply constructs the Fmutselneutral matrix for hyphy
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codon_dict  = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
codons      = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
purines     = ["A", "G"]
pyrims      = ["C", "T"]

def is_TI(source, target):
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False
        
def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff



def build_fmutsel_matrix(outfile):
    ''' Create FMutSel0 matrix.  ''' 
    #sij_set = set()
    matrix  = 'fmutsel = {61, 61, \n' # MG94
    
    for i in range(61):
        source = codons[i]
        sourceaa = codon_dict[source]
        for j in range(61):
            target = codons[j]
            targetaa = codon_dict[target]
            
            diff= get_nuc_diff(source, target)
            if len(diff) == 2:
                target_freq = "freq" + diff[1]   
#                 if sourceaa == targetaa:
#                     pr_fix = "1.0"
#                 else:
#                     pair = "".join(sorted(sourceaa + targetaa)).lower()
#                     if pair != (sourceaa+targetaa).lower():
#                         sij = "-1.0 * S" + pair
#                     else:
#                         sij = "S" + pair
#                         sij_set.add(sij)
#                     pr_fix = sij + "/(1-exp(-1*" + sij + "))"
# 
#                 # Create string for matrix element
#                 element = '{' + str(i) + ',' + str(j) + ',' + target_freq + "*" + pr_fix + "*t"  
                element = '{' + str(i) + ',' + str(j) + ',' + target_freq + "*t"  
                if is_TI(diff[0], diff[1]):
                    element += '*k'

                matrix += element + '}\n'
    # And save to file.
    with open(outfile, 'w') as outf:
        outf.write(matrix + '};\n')
    #for s in sij_set:
    #    print s
    #print len(sij_set)


build_fmutsel_matrix("fmutsel0.mdl")