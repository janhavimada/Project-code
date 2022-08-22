#Get acc id of refseq and position of cds. based on its stop codon being detected in which frame set the mumbers. seq_pos=end of cds, get align_pos and set frame.
import sys #Library for getting argument from command line
#Biopython for importing and manipulating sequence records
from Bio import SeqIO #For sequence records
from Bio.Seq import Seq #Manipulaing sequence
from Bio import Entrez #For fetching records
#Getting arguments from the command line
filename= sys.argv[1] #The name of the file which is the "species_segment" prefix
family=sys.argv[2] #Family name that the species belongs to
infile=f"../results/{family}/alignments/ali_{filename}.fasta" #Name of input file containing aligned fasta sequences
outfile=f"../results/{family}/orfs/ORF_{filename}.tsv" #Output file with all ORFs in TSV format
errorfile=f"../results/{family}/error/error_{filename}.txt" #File where the errors would be saved
references=sys.argv[3] #Path of the reference fasta file
Entrez.email = "ab@example.com"  # NCBI requires to set an email id

#Function that gets cds of required sequence with just its accession id
def findcds(acc_id):
    cds = {"+":[],"-":[]} #Dictionary with cds start and stop
    handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text") #Fetch genbank records for the reference
    for seq_record in SeqIO.parse(handle, "genbank"):
        for feature in (seq_record.features):
            if feature.type=="CDS":
                location = feature.location
                #CDS in the forward strand (+)
                if location.strand == 1:
                    cds["+"].append(f"{location.start+1}:{location.end}")
                #CDS in the reverse strand (-)
                if location.strand == -1:
                    cds["-"].append(f"{location.start+1}:{location.end}")
    return cds

#Function that finds sequence position with respect to alignment position or vice versa
def findpos(seq,start,fun):
    align_pos=0
    seq_pos=0
    for base in seq:
        #The alignment position (which incudes gaps) always increments by one
        align_pos += 1
        #If the base is not a gap then seq position increases
        if base != "-":
            seq_pos += 1
        #If seq position is entered then align_pos is returned
        if fun=="seq":
            if seq_pos==start:
                return align_pos
                break
        #If ali position is entered then seq_pos is returned
        if fun=="ali":
            if align_pos==start:
                return seq_pos
                break

#Function to return the orf in the given file
def getorf(sequence, seq_id, ali_start, ali_stop, file_out, strand, ref_check, cds_start, cds_stop):
    #If strand is -ve, the start and stop position have to be calculated again, according to the reverse complement
    if strand=="-":
        ali_start= len(sequence)-ali_start+1
        ali_stop= len(sequence)-ali_stop+1
    get_start=findpos(sequence, ali_start,"ali") #Start position of the cds according to the sequence
    #Presence of the bases before the cds starts to indicate presence or absence of UTR
    if get_start<=1:
        utr="N"
    else:
        utr="Y"
    set_frame= get_start%3 # The modulus helps determine frame of the cds start
    #The order of the frames is corrected according to the modulus, e.g if the cds starts at 2, then set_frame=2 and hence the frames would be -1(3), main(1), +1 (2) so on
    if set_frame==1:
        corrected_frame=[1,2,3]
    elif set_frame==2:
        corrected_frame=[3,1,2]
    elif set_frame==0:
        corrected_frame=[2,3,1]
    stop_codons = ["TAG", "TAA", "TGA"]
    newsequence = str(sequence).upper().replace("-","") #The original sequence without gaps
    for frame in range(3):
        seq_pos = frame #The first sequence would be 1, 2 or 3
        next_codon = "" #The variable stores the bases being read as triplets, like a codon
        for base in newsequence[frame:]:
            seq_pos += 1 #For each base read, sequence position goes up by 1
            next_codon += base
            #Once codon triplet is complete, the program will test if it is a stop codon, if a stop codon it determines the type of ORF
            if len(next_codon) == 3:
                if next_codon in stop_codons:
                    align_pos=findpos(sequence,seq_pos,"seq")
                    #If the position is the cds stop, or the frame is 1 and the position is greater than cds start, it is identified as an N-extension (Next)
                    if align_pos==ali_stop or (corrected_frame[frame]==1 and seq_pos>=get_start) :
                        type="Next"
                    #If the frame is 1 but stops before the cds starts that means it is interrupted within the 5' UTR, it is labelled "Inframe"
                    elif corrected_frame[frame]==1 and seq_pos<get_start:
                        type="Inframe"
                    #If neither then it has to be an overprinted ORF (UFO)
                    else:
                        type="UFO"
                    file_out.write(f"{corrected_frame[frame]}\t{strand}\t{type}\t{utr}\t{align_pos}\t{seq_pos}\t{seq_id}\t{next_codon}\t{newsequence[frame:seq_pos]}\t{ref_check}\t{cds_start}\t{cds_stop}\t{len(sequence)}\n")
                    #Exits the loop for current frame as a stop codon was found within it
                    break
                #If stop not found, resets the next_codon
                next_codon = ""

list_ids=[] #List containing accession ids of the sequences in the current species_segment input file
cds="" #Variable to contains the CDS dictionary
ref_check="NA" #The ref check is kept NA unless one of the sequences is in the reference fasta
#The sequences in the given fasta file are added to the list
for seq_record in SeqIO.parse(infile, "fasta"):
    seq_id=seq_record.id
    list_ids.append(seq_id)
#Checks the reference fasta file to see if any of its records is within the current species_segment file.
for ref_record in SeqIO.parse(references, "fasta"):
    ref_id= ref_record.id
    if ref_id in list_ids:
        ref_check=ref_id #Assigns it as the ref_check
        #Gets cds for the reference
        cds=findcds(ref_id)
        #Checks if the reference has a 5' UTR on both the strands, if not it skips this sequence and looks for the next reference
        try:
            if int(cds["+"][0].split(":")[0].replace("<",""))!=1:
                break
        except IndexError:
            pass
        try:
            if int(cds["-"][0].split(":")[0].replace(">",""))!=len(ref_record.seq):
                break
        except IndexError:
            pass
#If a reference sequence was not found it will choose the first sequence by default
if cds=="":
    cds = findcds(list_ids[0])
    ref_id= list_ids[0]
else:
    ref_id=ref_check #This is so that if the previous loop isnt terminated the last sequence shouldn't stored as ref_id

#At this point ref_id will store id of desired reference/first sequence
#ali_start and rev_ali_start are kept at 0 to check if the reference has appropriate cds
ali_start=0
rev_ali_start=0
#Will get the start and stop of the cds with respect to the alignment, so that it could be universally applied instead of fetching cds for every single sequences
for seq_record in SeqIO.parse(infile, "fasta"):
    #As cds is already known with respect to the reference sequence, once the reference sequence is extracted from the fasta file, the findpos function is used to get the cds start and stop with respect to the alignment
    if ref_id == seq_record.id:
        ref_seq= seq_record.seq
        rev_ref_seq= ref_seq.reverse_complement()
        #For positive strand
        try:
            cds["+"][0]
        except IndexError:
            pass
        else:
            try:
                cds_start=cds["+"][0].split(":")[0]
                ali_start=findpos(ref_seq,int(cds_start),"seq")
            except ValueError:
                ali_start=0
            cds_stop=cds["+"][0].split(":")[1]
            ali_stop=findpos(ref_seq,int(cds_stop.replace(">","")),"seq")
        #For negative strand
        try:
            cds["-"][0]
        except IndexError:
            pass
        else:
            try:
                rev_cds_start=cds["-"][0].split(":")[1]
                rev_ali_start=findpos(ref_seq,int(rev_cds_start),"seq")
            except ValueError:
                rev_ali_start=0
            rev_cds_stop=cds["-"][0].split(":")[0]
            rev_ali_stop=findpos(ref_seq,int(rev_cds_stop.replace(">","")),"seq")
        break
#If ali_start and rev_ali start is 0 that means either the sequences in the alignment are incomplete or none of them have a UTR, both making them unusable, hence the program is quit
if (ali_start+rev_ali_start)==0:
    with open (errorfile, "a") as file_out:
        file_out.write(f"The reference, {ref_id} has incomplete start sequence/ no UTRS")
    quit()
#Goes through each sequence in the fasta file, uses the getorf function to find and write the ORF
with open (outfile, "w") as file_out:
    file_out.write(f"Frame\tStrand\tType\tUTR\tAlignPos\tSeqPos\tSeqId\tCodon\tSeq\tRefSeq\tCdsStart\tCdsStop\tSeqLen\n") #Header to output file
    for seq_record in SeqIO.parse(infile, "fasta"):
        sequence_name= seq_record.description
        seq_id=seq_record.id
        sequence= seq_record.seq
        #For +ve strand
        try:
            cds["+"][0]
            if ali_start!=0:
                getorf(sequence, seq_id, ali_start, ali_stop, file_out, "+", ref_check, cds_start, cds_stop)
        except IndexError:
            pass
        #For -ve strand
        try:
            cds["-"][0]
            if rev_ali_start!=0:
                rev_sequence= sequence.reverse_complement()
                getorf(rev_sequence, seq_id, rev_ali_start, rev_ali_stop, file_out, "-", ref_check, rev_cds_start, rev_cds_stop)
        except IndexError:
            pass
