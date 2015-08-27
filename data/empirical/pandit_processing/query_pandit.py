# SJS
# Script to browse pandit and retrieve alignments.
# Pandit alignments are only considered if they meet criteria of: a) >=500 sequences, and b) <=0.5 average pairwise identity
# For those alignments, clean them up so that rows with >=25% data missing and columns with >=0.05% data are missing. If the final cleaning yields at least 450 sequences, save.

from Bio import AlignIO
from align_to_matrix import AlignMatrix
import subprocess
import urllib2
import re
import os

cull_columns    = 0.05 # Remove columns with >=5% of data missing, ambiguous
cull_rows       = 0.25 # Remove rows with >=25% data missing, ambiguous
pandit_numseq   = 500 # Only pull down PANDIT entries with >= 500 sequences
pandit_pairwise = 0.5 # Only pull down PANDIT entries with average pairwise identity <= 0.5
final_numseq    = 450 # Cleaned alignments are only saved if there are at least 450 sequences
final_numcol    = 150 # Cleaned alignments are only saved if there are at least 50 columns (since codon, len=150)


base_search_url = "http://www.ebi.ac.uk/goldman-srv/pandit/pandit.cgi?action=browse&key=REPLACE"
base_aln_url    = "http://www.ebi.ac.uk/goldman-srv/pandit/pandit.cgi?action=browse&fam=REPLACE&field=nam-dsq"
base_tree_url   = "http://www.ebi.ac.uk/goldman-srv/pandit/pandit.cgi?action=browse&fam=REPLACE&field=dtp-dph"

outdir = "pfam_queries/"


def clean_alignment(infile, outfile):
    
    a = AlignMatrix(alnfile = infile, format = "fasta", type = "codon")
    a.convert_aln_to_matrix()
    keep_cols = ( a.alnmat <= 0 ).sum(0) < cull_columns*a.numseq
    keep_cols[0] = True # always retain id
    a.alnmat = a.alnmat.compress(keep_cols, axis = 1)
    
    # Keep only rows with < 1% missing data or gaps
    keep_rows = ( a.alnmat[:,1:] <= 0 ).sum(1) < a.alnlen*cull_rows # gives the number of entries in each row which are bad
    a.alnmat = a.alnmat.compress(keep_rows, axis = 0)
    
    # Remove duplicates and save.
    a.convert_matrix_to_aln() 
    with open("rmdups.fasta", 'w') as outf:
        for k, v in a.aln_dict.items():
            outf.write(">" + k + "\n" + v + "\n")
    rmdup = subprocess.call("sh dup_hack.sh", shell=True)
    try:
        tempaln = AlignIO.read('rmdups.fasta.reduced', 'phylip-relaxed')
    except:
        tempaln = AlignIO.read('rmdups.fasta', 'fasta')
    
    # If final number of sequences is >= final_numseq, save it
    if len(tempaln) >= final_numseq and len(tempaln[0]) >= final_numcol:
        AlignIO.write(tempaln, outfile, "fasta")
        
    # Cleanup
    os.system("rm temp.fasta*")
    os.system("rm rmdups.fasta*")
    os.system("rm RAxML*")



def obtain_source(search_url, type):

    response = urllib2.urlopen(search_url)
    if type == "lines":
        page_source = response.readlines()
    elif type == "full":
        page_source = response.read()
        
    return page_source


def parse_pandit_page(query_pf):

    id_list = []
    
    # download the page source
    url = base_search_url.replace("REPLACE",query_pf)
    page_source = obtain_source(url, "lines")
    
    # parse the page source
    for line in page_source:
        if line.startswith('<tr bgcolor="#eeeeee">'):
            split_line = line.split('</td><td align="right">')[1:]
            numseq = int(split_line[0])
            pairid = float(split_line[2])
            
            # Grab the PFAM id if we are going to keep this 
            if numseq >= pandit_numseq and pairid <= pandit_pairwise:
                find_pfam = re.search("family\?entry=(PF\d+)", split_line[-1])
                if find_pfam:
                    pfam_id = find_pfam.group(1)
                    id_list.append(pfam_id)
    return id_list


    
def download_pandit_alignments(ids):
    
    for id in ids:
        #print "        Checking PFAM", id
        outfile = id + ".fasta"
        
        # download alignment page source
        url = base_aln_url.replace("REPLACE",id)
        page_source = obtain_source(url, "full")
        
        # Grab just the table lines
        page_source = page_source.split('<table cellspacing="0" cellpadding="4" width="100%" border="0">')[1].split('</table>')[0]
        
        # split table into key, value 
        page_source = page_source.replace("<tr>\n", "").replace("</tt></td>\n</tr>", "").replace('<td align="left" bgcolor="#eeeeee" nowrap><tt>', "").replace("</tt></td>", ";")
        page_source = page_source.replace("/", "_").replace(";\n", ";")
        pairs = page_source.split("\n")[1:-1]
        
        # Save to fasta file
        with open("temp.fasta", "w") as f:
            for entry in pairs:
                idseq = entry.split(";")
                seq = idseq[1].replace(".", "-").upper()
                f.write(">" + idseq[0] + "\n" + seq + "\n")
        
        print id
        # Clean and save if passes thresholds, and clean up
        clean_alignment("temp.fasta", outfile)


        
        
        
        
      
        
    
    
    
    
    
def main():

    for i in range(1,82):

        # Generate search key
        query_pf = str(i)
        while len(query_pf) != 3:
            query_pf = "0" + query_pf
        query_pf = "PF" + query_pf

        # download and parse the page source to obtain IDs we want
        ids = parse_pandit_page(query_pf)

        # grab the alignment and tree for each entry in ids
        download_pandit_alignments(ids)
    
main()
    
    

    

    