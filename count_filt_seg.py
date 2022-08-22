#Counting the number of segments that have passed all the filters and have been used for alignment and further final analysis
import sys #Library for getting argument from command line
import pandas as pd #Pandas for formating datatable
family=sys.argv[1] #Gets the family name that is being counted from the command line arguments
infile=f"../results/{family}/{family}_filtered.txt" #The txt file that contains the list of species segment and their counts for the family
outfile=f"../results/{family}/{family}_filt_segments.txt" #The file that will contain a more condensed tabular format of it
dataframe= pd.read_csv(infile,sep="\t") #Reads input file as TSV and stores as a data frame
pivoted = dataframe.pivot(index="Species", columns="Segment", values="Count") #Pivots the dataframe to condesne the format
pivoted.to_csv(outfile,sep="\t",na_rep=0) #Stores the condensed form to output file
