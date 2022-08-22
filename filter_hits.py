import sys #Library for getting argument from command line
#Biopython for importing and manipulating sequence records
from Bio import SeqIO
from Bio.Seq import Seq
#Accept and assign the arguments from command line
file=sys.argv[1]
family=sys.argv[2]
filename=f"../results/{family}/fastafiles/{file}.fasta" #Name of the fasta file containing the sequences
blast_txt=f"../results/{family}/blast_results/blast_{file}.txt" #Results from blast
seg=(file).split("_")[1].lower() #Segment name extracted from the filename
outfile=f"../results/{family}/error/error_{sys.argv[1]}.txt" #Output file that saves all the sequences that were chucked
sp=(file).split("_")[0].replace("-"," ")
rev=[] #List of id of the sequences that have to be reverse complemented
wrong=[] #List of sequences that do not match the segment hit
right=[] #List of sequences that match the segment hit
#Check blast hits to find the ones that align in opp direction, i.e. start>stop
with open (outfile,"w") as outwrite:
    with open (blast_txt,'r') as blast:
        for line in blast:
            line = line.split("\t")
            match=line[1] #The segment that it hit
            start=int(line[8]) #Get the start position in the orientation the hit matches
            stop=int(line[9]) #Get the stop position in the orientation the hit matches
            acc_id=line[0] #Get the accession id of the sequence
            #If start> stop that means that the hit aligns in reverse orietation of the sequence, hence the sequence has to be reverse complemented
            if start>stop:
                rev.append(acc_id)
            #If the segment name matches the hit then the accession id is added to the "right" list
            if match==seg:
                right.append(acc_id)
            #If the segment name doesn't match the hit then accession id is added to the "wrong" list with the name of the segment that it hits
            if match!=seg:
                wrong.append(acc_id)
                outwrite.write(acc_id+"\t"+match+"\n")
    record=[] #The updated fasta records will be added to this list
    #The fasta file and the error file is updated
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id in right:
            #The sequences that are right will be added to the fasta file and will be reverse complemented before adding if they are in rev
            if seq_record.id in rev:
                my_seq=  (seq_record.seq).reverse_complement()
                seq_record.seq= my_seq
            record.append(seq_record)
        #If the record is not in "right" or "wrong" it might not have any hits, so it is added to the error file
        if seq_record.id not in right:
            if seq_record.id not in wrong:
                outwrite.write(seq_record.id+"\t no hit \n")
    SeqIO.write(record,filename,"fasta")
#Writes the count of the filtered sequences to be used for analysis into the text file
with open(f"../results/{family}/{family}_filtered.txt","a") as segfile:
    segfile.write(f"{sp}\t{seg}\t{len(right)}\n")
