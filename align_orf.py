import sys
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
tsvin=sys.argv[1]
# orf table has the sequences, get the sequences and translate and align.
# create tsv with only interesting orfs.
# get the name of the sequence, the align pos and frame,
# open orf file: get the sequences, ref seq should be on top
with open (tsvin,"r") as tsv:
    for line in tsv:
        splits=line.split("\t")
        family=splits[0]
        species=splits[1].replace(" ","-")
        segment=splits[2]
        strand=splits[3]
        frame=int(splits[4])
        type=splits[5]
        utr=splits[6]
        alnstop=int(splits[7])
        count=int(splits[8])
        total=int(splits[9])
        proportion=splits[10]
        ref_seq=splits[11]
        cds_start=splits[12]
        cds_stop= splits[13]
        seq_length=int(splits[14])
        stop_pos=int(splits[15])
        Nseq=splits[16]
        AAseq=splits[17]
        AAlen=int(splits[18].replace("\n",""))
        seq_dict={ref_seq:""}
        with open(f"../results/{family}/orfs/ORF_{species}_{segment}.tsv") as orffile:
            in_file=f"../results/prot/{species}_{segment}_{frame}.fasta"
            counter=0
            #create a fasta file which can be input file for clustal omega
            for line in orffile:
                if counter==0:
                    counter+=1
                    continue
                splits=line.split("\t")
                seq_id=splits[6]
                if int(splits[4])==alnstop and int(splits[0])==frame and splits[1]==strand:
                    AAseq=Seq(splits[8][:-3]).translate()
                    record= SeqRecord(AAseq,id=seq_id)
                    seq_dict[seq_id]= record
            seq_records= list(seq_dict.values())
            SeqIO.write(seq_records, in_file, "fasta")
            out_file=f"../results/prot/ali_{species}_{segment}_{frame}.fasta"
            clus_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,verbose=True, auto=True, force=True)
            clus_cline()
