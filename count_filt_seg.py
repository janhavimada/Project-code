#Counting the number of segments that have passed all the filters and have been used for alignment and further final analysis
import sys #Library for getting argument from command line
#Biopython for importing and manipulating sequence records
import pandas as pd
#Accept and assign the arguments from command line
family=sys.argv[1] #Family name that is being counted
infile=f"{family}/{family}_filtered.txt" #The txt file that contains the list of species segment and their counts for the family
outfile=f"{family}/{family}_filt_segments.txt" #The file that will contain a more condensed tabular format of it

dataframe= pd.read_csv(infile,sep="\t")
pivoted = dataframe.pivot(index="Species", columns="Segment", values="Count")
pivoted.to_csv(outfile,sep="\t",na_rep=0)
