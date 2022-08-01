# Predicting UFOs and N-ext by cap snatching
The following code allows user to give segmented negative strand RNA virus sequences and obtain possible UFOs and N-ext of cannonical proteins.

## Input
The scripts should be downloaded and saved in a sub-directory called "code" within the directory that the results are expected in. The input data is also stored in a separate sub-directory called "data".
Eg:
!(/directories-img.png)
Three types of input data are expected from the user.
1. CSV containing metadata for all the sequences, including RefSeq sequences. Ideally downloaded from the NCBI Virus database. It should be in the format of:
2. Fasta file of the all the sequences from above said metadata.
3. Another fasta file that only contains the sequences that have RefSeq annotation. These should also contain the segment name at the end of description

## Usage
The `proj_env.yml` contains the environment. It can be activated by the following code:
```
conda env create -n proj_env --file proj_env.yml
conda activate proj_env
```
The `telescope.sh` is the executable. The other scripts do not have to be run separately but need to be in the same directory as the executable.
```
./telescope.sh $1 $2 $3
```
Here $1, $2 and $3 are placeholders for the input given as mentioned in above section. Using the exampke data the code would look like:
```
./telescope.sh ../data/polyploviricotina_seq.csv ../data/polyplo_sequences.fasta ../data/ref_seqog.fasta
```
The original data is lengthy and would take an hour to run.
All the results will appear in the "results" subdirectory within the same directory as the code. The results appear separately for all the families.
"ORFs_mastertable.tsv" is useful to get a summarised version.
