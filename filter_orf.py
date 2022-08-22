#Filter the orfs by removing shorter sequences, remove redundancy and get a list of proteins
import sys
#Biopython for importing and manipulating sequence records
from Bio import SeqIO #For sequence records
from Bio.Seq import Seq #Manipulaing sequence
filename= sys.argv[1] #The name of the file which is the "species_segment" prefix
family=sys.argv[2] #Family name that the species belongs to
reffasta= sys.argv[3] #Path of the reference fasta file
species=filename.split("_")[0].replace("-"," ") #Gets species name from  the filename
segment=filename.split("_")[1] #Gets segment name from  the filename
infile=f"../results/{family}/orfs/ORF_{filename}.tsv" #Input file with all ORFs in TSV format
count_file=f"../results/{family}/{family}_filtered.txt" #Text file containg the total counts analysed
outfile="../results/ORFs_mastertable.tsv" #Main output file containing the ORFs for all the species and segments in a condensed format
#A nested dictionary containing the different ORFs based on stop positions by alignment as per the Strand > Type > Presence/absence of UTR
stop_pos={"+":{"UFO":{"Y":{},"N":{}},"Next":{"Y":{},"N":{}},"Inframe":{"Y":{},"N":{}}},"-":{"UFO":{"Y":{},"N":{}},"Next":{"Y":{},"N":{}},"Inframe":{"Y":{},"N":{}}}}
#Creates list of the IDs of sequences in the reference fasta file
list_of_ref=[]
for ref_record in SeqIO.parse(reffasta, "fasta"):
    list_of_ref.append(ref_record.id)
#Opens input file to get the ORFs and add them to the disctionary
with open(infile,"r") as orftsv:
    counter=0
    for line in orftsv:
        if counter==0:
            counter+=1
            continue
        splits = line.split("\t")
        alnstop= int(splits[4]) #The stop position by alignment
        frame= int(splits[0]) #The frame of ORF
        strand= splits[1] #The strand on which ORF is present: +/- ve
        type= splits[2] #The type of ORF: UFO/Next/Inframe
        utr= splits[3] #If UTR is absent/present: N/Y
        seqpos= int(splits[5]) #The stop position as per sequence
        seqid= splits[6] #The accession id of the sequence
        Nseq= splits[8] #The nucleotide sequence of ORF
        ref_seq= splits[9] #Contents of reference column
        cds_start=splits[10] #CDS start position as per reference/first sequence
        cds_stop=splits[11] #CDS stop position as per reference/first sequence
        seq_len=splits[12].replace("\n","") #Length of the nucleotide sequence
        #Checks if the record for the given ORF exists, if not adds the record
        #Record is in the format of list with [frame, count, reference, accession, cds start, cds stop, length, List of accesion IDS,  a list for first sequence-[n seq, aa seq, length, position], a list for reference seq-[n seq, aa seq, length, position],]
        if alnstop not in stop_pos[strand][type][utr]:
            AAseq= str(Seq(Nseq).translate()).replace("*","") #Translates nucleotide sequence
            stop_pos[strand][type][utr][alnstop]=[frame, 0, ref_seq, seqid, cds_start, cds_stop, seq_len, [], [Nseq, AAseq, len(AAseq),seqpos]]
        if alnstop in stop_pos[strand][type][utr]:
            stop_pos[strand][type][utr][alnstop][1]+=1 #Ups the count by one everytime it comes across this ORF
            stop_pos[strand][type][utr][alnstop][7].append(seqid) #Adds accesion ID for each sequence having the ORF
        if seqid==stop_pos[strand][type][utr][alnstop][2]:
            AAseq= str(Seq(Nseq).translate()).replace("*","") #Translates nucleotide sequence
            stop_pos[strand][type][utr][alnstop].append([Nseq, AAseq, len(AAseq),seqpos]) #Adds details of reference to end of list
#To get the count of number of sequence analysed
with open(count_file,"r") as readfile:
    for line in readfile:
        sp=line.split("\t")[0]
        seg=line.split("\t")[1]
        if sp==species and seg==segment:
            filt_count=int(line.split("\t")[2])
            break
#Writes from the dictionary to the output master TSV
with open(outfile,"a") as outwrite:
    for strand in stop_pos:
        for type in stop_pos[strand]:
            for utr in stop_pos[strand][type]:
                for alnstop in stop_pos[strand][type][utr]:
                    dict=stop_pos[strand][type][utr][alnstop]
                    frame=dict[0]
                    count=dict[1]
                    ref_seq=dict[2]
                    cds_start=dict[4]
                    cds_stop= dict[5]
                    seq_len= dict[6]
                    acc_id= ",".join(dict[7])
                    #Checks if the ORF is also seen in the sequence used as reference
                    try:
                        ref_list=dict[9]
                    #If not the first sequence with the ORF is used as reference
                    except IndexError:
                        ref_list=dict[8]
                        ref_seq=dict[3]
                    #Checks if the sequence used as the reference is actually in the reference fasta file
                    if ref_seq in list_of_ref:
                        is_ref="Y"
                    else:
                        is_ref="N"
                    Nseq=ref_list[0]
                    AAseq=ref_list[1]
                    AAlen=ref_list[2]
                    seqpos=ref_list[3]
                    outwrite.write(f"{family}\t{species}\t{segment}\t{strand}\t{frame}\t{type}\t{utr}\t{alnstop}\t{count}\t{filt_count}\t{count/filt_count}\t{ref_seq}\t{is_ref}\t{cds_start}\t{cds_stop}\t{seq_len}\t{seqpos}\t{Nseq}\t{AAseq}\t{AAlen}\t{acc_id}\n")
