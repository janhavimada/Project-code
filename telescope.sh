#! /bin/bash
metadata=$1 #"../data/polyploviricotina_seq.csv"
fastafile=$2 #"../data/polyplo_sequences.fasta"
ogreffasta=$3 #"../data/ref_seqog.fasta"
reffasta="../data/ref_seq.fasta"
mkdir -p ../results
mkdir -p ../results/prot
#Create reference database
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
    #Perform BLAST to cleanse the data. i.e. segment should match reference.
    blastn -db REF -max_target_seqs 1 -max_hsps 1 -query ../results/${family}/fastafiles/${file}.fasta -out ../results/${family}/blast_results/blast_${file}.txt -outfmt 6 2>/dev/null
    #Filter the blast results
    python3 filter_hits.py ${file} ${family}
    #Remove empty fastafiles and skip the ssegment if empty
    if ! [ -s ../results/${family}/fastafiles/${file}.fasta ]; then rm -f ../results/${family}/fastafiles/${file}.fasta; continue; fi
    #Perform MAFFT to align the remaining sequences in the fasta file
    mafft --quiet ../results/${family}/fastafiles/${file}.fasta> ../results/${family}/alignments/ali_${file}.fasta
    #Use script to detect the position of stop codons
    python3 getstop.py ${file} ${family} ${ogreffasta}
    #Remove empty error files
    if ! [ -s ../results/${family}/error/error_${file}.txt ]; then rm -f ../results/${family}/error/error_${file}.txt; fi
    #Remove empty ORF result files and skip the segment if empty
    if ! [ -s ../results/${family}/orfs/ORF_${file}.tsv ]; then rm -f ../results/${family}/orfs/ORF_${file}.tsv; continue; fi
    #Add a summary to the Master table
    python3 filter_orf.py ${file} ${family}
  done
  #Convert the structure of the filtered segments table using pandas
  python3 count_filt_seg.py ${family}
 done
