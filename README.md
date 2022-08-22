# Predicting UFOs and N-ext by cap snatching
The following code allows user to give segmented negative strand RNA virus sequences and obtain possible UFOs and N-ext of cannonical proteins.

## Dependencies
The `proj_env.yml` contains the environment. It can be activated by the following code:
```
conda env create -n proj_env --file proj_env.yml
conda activate proj_env
```
The command line programs that must be installed for the pipeline to run are as follows:
MAFFT v7.505 https://mafft.cbrc.jp/alignment/software/
BLAST+ v2.13.0 https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

If `align_orf.py` (not a part of the pipeline, can be used for preparing amino acid alignment files for protein structure prediction) needs to be run then Clustal Omega v1.2.3 is required. http://www.clustal.org/omega/#Download

## Input
The scripts should be downloaded and saved in a sub-directory called "code" within the directory that the results are expected in. The input data is also stored in a separate sub-directory called "data".Eg:

![](/directories-img.png)

Three types of input data are expected from the user.
1. CSV containing metadata for all the sequences, including RefSeq sequences. Ideally downloaded from the NCBI Virus database. It should be in the format of:
Accession, Species, Genus, Family, Molecule_type (DNA/RNA), Length, Segment, GenBank_Title
2. Fasta file of the all the sequences from above said metadata.
3. Another fasta file that only contains the sequences that are being used as reference. These should also contain the segment name at the end of description.
Example dataset is available in `data.zip`. `polyploviricotina_seq.csv` `polyplo_sequences.fasta` `ref_seqog.fasta` are the input files 1,2,3 respectively.
## Usage
The `telescope.sh` is the executable. The other scripts do not have to be run separately but need to be in the same directory as the executable.
```
./telescope.sh $1 $2 $3
```
Here $1, $2 and $3 are placeholders for the input given as mentioned in above section. For the example dataset it would be:
```
./telescope.sh ../data/polyploviricotina_seq.csv ../data/polyplo_sequences.fasta ../data/ref_seqog.fasta
```
Depending on the amount of data it might take up to an hour to run.
All the results will appear in the "results" subdirectory within the same directory as the code. The results appear separately for all the families.
"ORFs_mastertable.tsv" is useful to get a summarised version.
