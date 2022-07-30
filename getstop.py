#Get acc id of refseq and position of cds. based on its stop codon being detected in which frame set the mumbers. seq_pos=end of cds, get align_pos and set frame.
import sys #Library for getting argument from command line
from Bio import SeqIO #For sequence records
from Bio.Seq import Seq
from Bio import Entrez #For fetching records
#Getting arguments from the command line
filename= sys.argv[1]
family=sys.argv[2]
infile=f"../results/{family}/alignments/ali_{filename}.fasta"
outfile=f"../results/{family}/orfs/ORF_{filename}.tsv"
errorfile=f"../results/{family}/error/error_{filename}.txt"
references=sys.argv[3]
Entrez.email = "ab@example.com"  # NCBI requires to set an email id

#Function that gets cds of required seq with its accession id
def findcds(acc_id):
    cds = {"+":[],"-":[]} #Dictionary with cds start and stop
    handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text") #Fetch records for the ref_seq
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

#Function that finds seq position with respect to align position or vice versa
def findpos(seq,start,fun):
    align_pos=0
    seq_pos=0
    for base in seq:
        # the alignment position (which incudes gaps) always increments by one
        align_pos += 1
        # if the base is not a gap then seq position increases
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
    if strand=="-":
        ali_start= len(sequence)-ali_start+1
        ali_stop= len(sequence)-ali_stop+1
    get_start=findpos(sequence, ali_start,"ali")
    set_frame= get_start%3
    if get_start<=1:
        utr="N"
    else:
        utr="Y"
    # get_stop= findpos(sequence, ali_stop,"ali")
    if set_frame==1:
        corrected_frame=[1,2,3]
    elif set_frame==2:
        corrected_frame=[3,1,2]
    elif set_frame==0:
        corrected_frame=[2,3,1]
    stop_codons = ["TAG", "TAA", "TGA"]
    newsequence = sequence.upper().replace("-","")
    for frame in range(3):
        seq_pos= frame
        next_codon=""
        for base in newsequence[frame:]:
            seq_pos += 1
            next_codon += base
            # print(next_codon + " " + str(align_pos) + " " + str(seq_pos))
            if len(next_codon) == 3:
                if next_codon in stop_codons:
                    align_pos=findpos(sequence,seq_pos,"seq")
                    if align_pos==ali_stop or (corrected_frame[frame]==1 and seq_pos>=get_start) :
                        type="Next"
                    elif corrected_frame[frame]==1 and seq_pos<=get_start:
                        type="Inframe"
                    else:
                        type="UFO"
                    file_out.write(f"{corrected_frame[frame]}\t{strand}\t{type}\t{utr}\t{align_pos}\t{seq_pos}\t{seq_id}\t{next_codon}\t{newsequence[frame:seq_pos]}\t{ref_check}\t{cds_start}\t{cds_stop}\t{len(sequence)}\n")
                    # exit the base loop as found a codon so no need to look further
                    break
                    # reset the next_codon as hit length three
                next_codon = ""

# USE BIOPYTHON TO READ IN A FASTA FILE (CAN BE ALIGNED OR UNALIGNED DOESN"T MATTER)
# LOOP - ONE YOU HAVE A SEQUENCE YOU CAN RUN THE BELOW ON EACH
# YOU CAN THEN REVERSE COMP THE SEQUENCE AND RE_RUN AS WELL
list_ids=[]
cds=""
ref_check="NA"
for seq_record in SeqIO.parse(infile, "fasta"):
    seq_id=seq_record.id
    list_ids.append(seq_id)
for ref_record in SeqIO.parse(references, "fasta"):
    ref_id= ref_record.id
    if ref_id in list_ids:
        ref_check=ref_id
        cds=findcds(ref_id)
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

if cds =="":
    cds = findcds(list_ids[0])
    ref_id= list_ids[0]
        #Find cds save the start, stop in a variable.
        #If -, reverse comp
        #Find align pos. thus frame.
#At this point ref_id will store id of desired ref seq
ali_start=0
rev_ali_start=0
for seq_record in SeqIO.parse(infile, "fasta"):
    if ref_id == seq_record.id:
        ref_seq= seq_record.seq
        rev_ref_seq= ref_seq.reverse_complement()
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

if (ali_start+rev_ali_start)==0:
    with open (errorfile, "a") as file_out:
        file_out.write(f"The reference, {ref_id} has incomplete start sequence/ no UTRS")
    quit()

with open (outfile, "w") as file_out:
    file_out.write(f"Frame\tStrand\tType\tUTR\tAlignPos\tSeqPos\tSeqId\tCodon\tSeq\tRefSeq\tCdsStart\tCdsStop\tSeqLen\n")
    for seq_record in SeqIO.parse(infile, "fasta"):
        sequence_name= seq_record.description
        seq_id=seq_record.id
        sequence= seq_record.seq
        try:
            cds["+"][0]
            if ali_start!=0:
                getorf(sequence, seq_id, ali_start, ali_stop, file_out, "+", ref_check, cds_start, cds_stop)
        except IndexError:
            pass
        try:
            cds["-"][0]
            if rev_ali_start!=0:
                rev_sequence= sequence.reverse_complement()
                getorf(rev_sequence, seq_id, rev_ali_start, rev_ali_stop, file_out, "-", ref_check, rev_cds_start, rev_cds_stop)
        except IndexError:
            pass
# print(Seq(ref_seq[914:1658]).reverse_complement())
