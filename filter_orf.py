#Filter the orfs by removing shorter sequences, remove redundancy and get a list of proteins
import sys
from Bio import SeqIO
from Bio.Seq import Seq
#CSV will be like, SPECIES | SEGMENT | ORF LENGTH | AA LENGTH | EXT TYPE (n ext or diff frame)|
filename= sys.argv[1]
family=sys.argv[2]
reffasta= sys.argv[3]
species=filename.split("_")[0].replace("-"," ")
segment=filename.split("_")[1]
infile=f"../results/{family}/orfs/ORF_{filename}.tsv"
count_file=f"../results/{family}/{family}_filtered.txt"
outfile="../results/ORFs_mastertable.tsv"
stop_pos={"+":{"UFO":{"Y":{},"N":{}},"Next":{"Y":{},"N":{}},"Inframe":{"Y":{},"N":{}}},"-":{"UFO":{"Y":{},"N":{}},"Next":{"Y":{},"N":{}},"Inframe":{"Y":{},"N":{}}}}
list_of_ref=[]
for ref_record in SeqIO.parse(reffasta, "fasta"):
    list_of_ref.append(ref_record.id)
with open(infile,"r") as orftsv:
    counter=0
    for line in orftsv:
        if counter==0:
            counter+=1
            continue
        splits = line.split("\t")
        alnstop= int(splits[4])
        frame= int(splits[0])
        strand= splits[1]
        type= splits[2]
        utr= splits[3]
        seqpos= int(splits[5])
        seqid= splits[6]
        Nseq= splits[8]
        ref_seq= splits[9]
        cds_start=splits[10]
        cds_stop=splits[11]
        seq_len=splits[12].replace("\n","")
        if alnstop not in stop_pos[strand][type][utr]:
            AAseq= str(Seq(Nseq).translate()).replace("*","")
            stop_pos[strand][type][utr][alnstop]=[frame, 0, ref_seq, seqid, cds_start, cds_stop, seq_len, [], [Nseq, AAseq, len(AAseq),seqpos]]
        if alnstop in stop_pos[strand][type][utr]:
            stop_pos[strand][type][utr][alnstop][1]+=1
            stop_pos[strand][type][utr][alnstop][7].append(seqid)
        if seqid==stop_pos[strand][type][utr][alnstop][2]:
            AAseq= str(Seq(Nseq).translate()).replace("*","")
            stop_pos[strand][type][utr][alnstop].append([Nseq, AAseq, len(AAseq),seqpos])
with open(count_file,"r") as readfile:
    for line in readfile:
        sp=line.split("\t")[0]
        seg=line.split("\t")[1]
        if sp==species and seg==segment:
            filt_count=int(line.split("\t")[2])
            break

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
                    try:
                        ref_list=dict[9]
                    except IndexError:
                        ref_list=dict[8]
                        ref_seq=dict[3]
                    if ref_seq in list_of_ref:
                        is_ref="Y"
                    else:
                        is_ref="N"
                    Nseq=ref_list[0]
                    AAseq=ref_list[1]
                    AAlen=ref_list[2]
                    seqpos=ref_list[3]
                    outwrite.write(f"{family}\t{species}\t{segment}\t{strand}\t{frame}\t{type}\t{utr}\t{alnstop}\t{count}\t{filt_count}\t{count/filt_count}\t{ref_seq}\t{is_ref}\t{cds_start}\t{cds_stop}\t{seq_len}\t{seqpos}\t{Nseq}\t{AAseq}\t{AAlen}\t{acc_id}\n")
