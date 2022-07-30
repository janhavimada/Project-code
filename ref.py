#Code for creating reference fasta file for the database. It includes RefSeq sequences
import sys #Library for getting argument from command line
#Biopython for importing and manipulating sequence records
from Bio import SeqIO
from Bio.Seq import Seq
infile=sys.argv[1] #Filename is accepted from the command line
outfile=sys.argv[2]
record=[] #The updated fasta records
#The segment name has to appear before the accession id so that BLAST will be easier
for seq_record in SeqIO.parse(infile, "fasta"):
    #Extract segment name, convert to lower case, and remove spaces, hyphens, underscores for consistency and avoiding duplicates
    seg=(((seq_record.description.split("|")[2].replace(" ","")).replace("-","")).replace("_","")).lower()
    #Some segments have been inconcistenly named
    if seg=="s(small)":
        seg="s"
    if seg=="small":
        seg="s"
    if seg=="large":
        seg="l"
    if seg=="medium" or seg=="middle":
        seg="m"
    seq_record.id=seg+ " "+seq_record.id
    record.append(seq_record) #Updated record saved to the list
SeqIO.write(record,outfile,"fasta") #Updated records written into fastafile
