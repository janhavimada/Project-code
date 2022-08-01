#! /bin/bash
metadata=$1 #"../data/polyploviricotina_seq.csv"
fastafile=$2 #"../data/polyplo_sequences.fasta"
ogreffasta=$3 #"../data/ref_seqog.fasta"
reffasta="../data/ref_seq.fasta"
#create REF
mkdir -p ../results
echo "Creating reference database..."
python3 ref.py ${ogreffasta} ${reffasta} #Code for creating reference fasta file
makeblastdb -in ${reffasta} -dbtype nucl -out REF
echo "Created ref db"
#Create list of family:
tail -n +2 ${metadata} | awk -F',' '{print $4}' | sort -u > ../data/familynames.txt
#Creates file for final output with headers
echo  $'Family\tSpecies\tSegment\tStrand\tFrame\tType\tUTR\tCount\tTotal\tProportion\tRefSeq\tCDS Start\tCDS Stop\tSeq Length\tStop Pos\tN seq\tAA seq\tAA len' > ../results/ORFs_mastertable.tsv
for family in $(echo $(cat ../data/familynames.txt))
do
  echo "Processing ${family} ..."
  mkdir -p ../results/${family}
  mkdir -p ../results/${family}/fastafiles
  mkdir -p ../results/${family}/blast_results
  mkdir -p ../results/${family}/error
  mkdir -p ../results/${family}/alignments
  mkdir -p ../results/${family}/orfs
  #Process metadata to separate csvs by family
	grep ${family} ${metadata} > ../results/${family}/${family}.csv
  #Get acc_id, species and thus seperate fasta files using the python script
  echo "Creating fasta files for ${family}"
	python3 count_seg.py ${family} ${fastafile} #Counting the segments for each species and making summary, splitting the fasta files
  #Loop that goes over the fasta files in family.txt
  for file in $(echo $(cat ../results/${family}/${family}.txt))
  do
    echo "processing ${file}"
    blastn -db REF -max_target_seqs 1 -max_hsps 1 -query ../results/${family}/fastafiles/${file}.fasta -out ../results/${family}/blast_results/blast_${file}.txt -outfmt 6 2>/dev/null  #BLAST
    python3 filter_hits.py ${file} ${family}
    if ! [ -s ../results/${family}/fastafiles/${file}.fasta ]; then rm -f ../results/${family}/fastafiles/${file}.fasta; continue; fi
    # echo "filtered blast for ${file}"
    mafft --quiet ../results/${family}/fastafiles/${file}.fasta> ../results/${family}/alignments/ali_${file}.fasta
    # echo "aligned ${file}"
    python3 getstop.py ${file} ${family} ${ogreffasta}
    # echo "got orf for ${file}"
    if ! [ -s ../results/${family}/error/error_${file}.txt ]; then rm -f ../results/${family}/error/error_${file}.txt; fi
    if ! [ -s ../results/${family}/orfs/ORF_${file}.tsv ]; then rm -f ../results/${family}/orfs/ORF_${file}.tsv; continue; fi
    python3 filter_orf.py ${file} ${family}
  done
 done
