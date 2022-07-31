# Predicting UFOs and N-ext by cap snatching
The following

## Input
The scripts should be downloaded and saved in a sub-directory within the directory that the results are expected in. The input data is also stored in a separate sub-directory called "data". Eg:

Three types of input data are expected from the user.
1. CSV containing metadata for all the sequences, including RefSeq sequences. Ideally downloaded from the NCBI Virus database. It should be in the format of:

2. Fasta file of the all the sequences from above said metadata.

3. Another fasta file that only contains the sequences that have RefSeq annotation. These should also contain the segment name at the end of description

## Usage
The ./cumulative.sh is the executable. The other scripts do not  have to be run separately.
All the results will appear in the "results" subdirectory which is in the same directory as the code. The results appear separately for all the families.
"ORFs_mastertable.tsv" is a summarised results table.
