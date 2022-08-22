#Optional code, not part of the pipeline
#Prepares aligned fasta files of the amino acid sequences from the ORFs which can be used in ColabFold for sequence prediction
#However, has to be converted into A3M for ColabFold using FormatSeq
import sys #Library for getting argument from command line
#Biopython for importing and manipulating sequence records
from Bio.Align.Applications import ClustalOmegaCommandline #For MSA of amino acids
from Bio.SeqRecord import SeqRecord #Create sequence records for the fasta file
from Bio.Seq import Seq #To translate nucleotide sequnces
from Bio import SeqIO #Read and write fastafiles
tsvin=sys.argv[1] #TSV with the ORFs. This can be created by selecting only desired sequences from the master table, without the header
#From TSV, get the name of the family, species, segment, frame, the stop position by alignment and the reference sequence
with open (tsvin,"r") as tsv:
    for line in tsv:
        splits=line.split("\t")
        family=splits[0]
        species=splits[1].replace(" ","-")
        segment=splits[2]
        strand=splits[3]
        frame=int(splits[4])
        alnstop=int(splits[7])
        ref_seq=splits[11]
        #Dictionary of the sequences, uses the reference as the first sequence so that it appears first on the fasta file
        if ref_seq !="NA":
            seq_dict={ref_seq:""}
        else:
            seq_dict={}
        #Opens orf tsv file to get the sequences having the ORF
        with open(f"../results/{family}/orfs/ORF_{species}_{segment}.tsv") as orffile:
            in_file=f"../results/prot/{species}_{segment}_{frame}.fasta" #FIle that stores all sequences in fasta format for Clustal Omega
            counter=0
            for line in orffile:
                #Skips header
                if counter==0:
                    counter+=1
                    continue
                splits=line.split("\t")
                seq_id=splits[6]
                #Checks if the ORF is the same as the one desired
                if int(splits[4])==alnstop and int(splits[0])==frame and splits[1]==strand:
                    AAseq=Seq(splits[8][:-3]).translate() #Translates the amino acid sequence from the ORF
                    record= SeqRecord(AAseq,id=seq_id) #Stores the AA sequence and ID as per biopython SeqRecord object format
                    seq_dict[seq_id]= record #Adds record to dictionary
            seq_records= list(seq_dict.values()) #Converts the sequence dictionary to a list
            SeqIO.write(seq_records, in_file, "fasta") #Writes amino acid sequences into the fasta file
            out_file=f"../results/prot/ali_{species}_{segment}_{frame}.fasta" #Output fasta file from Clustal Omega containg aligned sequences
            #Command Line interface for ClustalOmega
            clus_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,verbose=True, auto=True, force=True)
            clus_cline()
